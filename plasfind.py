import argparse, io, shutil, os
import subprocess as sp
import pandas as pd
from Bio import SeqIO
from codecs import decode
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram


import matplotlib.pyplot as plt
import umap.umap_ as umap

def arguments():
    parser = argparse.ArgumentParser(description="Finds Pantoea plasmids based on homology")
    parser.add_argument('-d', '--plasdb', default="PlasDB/Plasmids", type=str, help="DB containing plasmid fasta")
    parser.add_argument('-e', '--plasdb_ext', default=".fasta", type=str, help="extension to access the sequences in plasmid database")
    parser.add_argument('-i', '--infile', required=True, type=str, help="Fasta with contigs to be classified fasta")
    parser.add_argument('-c', '--chromdb', default="chromdb/chromosome", type=str, help="DB containing plasmid fasta")
    parser.add_argument('-t', '--tempdir', default="tmp", type=str, help="Directory to contain temporary files.")
    parser.add_argument('-n', '--num_threads', default=32, type=int, help="Number of threads for blast searches.")
    parser.add_argument('-o', '--outputdir', default="output", type=str, help="output directory")


    args = parser.parse_args()
    return args


def blast_run(query, db_path, perc=80, eval=0.00005, threads=1):
    cmd=f"blastn -query {query} -perc_identity {perc} -evalue {eval} -num_threads {threads} " \
        f"-outfmt '10 sseqid qlen slen length qcovus pident evalue qstart qend' -db {db_path}"
    proc1=sp.run(cmd, shell=True, stdout=sp.PIPE)
    return proc1.stdout

def blast2df(result):
    col_names=["qtitle", "query_length", "subject_length", "alignment_length", "query_coverage", "percent_identical", "e value", "qstart", "qend"]
    df=pd.read_csv(
        io.StringIO(decode(result)),
        names=col_names
    )
    df["Match_base_pid"]=df["query_coverage"]*df["percent_identical"]/100
    df=df.drop_duplicates(subset="qtitle")
    return df

def read_fasta(fastafile):
    '''reads fasta file'''
    recordlist = []
    #open file
    with open(fastafile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            recordlist.append(record)
    return recordlist

def read_fasta_single(fastafile):
    with open(fastafile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return record

def write_fasta(path, recordlist):
    with open(path, "w") as output_handle:
        SeqIO.write(recordlist, output_handle, "fasta")

def makedir(path: str, force_make=False):
    try:
        os.mkdir(path)
        return path
    except:
        if not force_make:
            return path
    i = 1
    if force_make:
        while True:
            new_name =f"{path}({i})"
            try:
                os.mkdir(new_name)
                break
            except:
                i+=1
        return new_name

def fasta_breaker(fasta, tempdir):
    records=read_fasta(fasta)
    i=1
    for contig in records:
        fasta_location=f"{tempdir}/contig_{i}"
        write_fasta(fasta_location, contig)
        i+=1

def contig_class(contigs_loc, blast_plas_loc, blast_chrom_loc, plasdb, chromdb, threads):
    contig_data_arr=[]
    for contig in os.listdir(contigs_loc):
        contig_path=f"{contigs_loc}/{contig}"
        blast_plas_path=f"{blast_plas_loc}/{contig}.csv"
        blast_chrom_path=f"{blast_chrom_loc}/{contig}.csv"

        plas_df=blast2df(blast_run(contig_path, plasdb, threads=threads))
        plas_df.to_csv(blast_plas_path)

        chrom_df=blast2df(blast_run(contig_path, chromdb, threads=threads))
        chrom_df.to_csv(blast_chrom_path)

        plas_result=isplasmid(plas_df, chrom_df)

        contig_data_arr.append((contig_path, chrom_df, plas_df, plas_result))
    return contig_data_arr

def write_chrom_fasta(contig_arr, location):
    chromosome_records=[]
    for contig in contig_arr:
        if not contig[3]:
            record=read_fasta_single(contig[0])
            chromosome_records.append(record)
    SeqIO.write(chromosome_records, location, "fasta")


def write_txt(path, text):
    with open(path, "w") as f:
        f.write(text)

def isplasmid(plas_df, chrom_df):
    max_plas = plas_df["query_coverage"].max()
    if len(plas_df) == 0:
        max_plas = 0


    max_chrom = chrom_df["query_coverage"].max()
    if len(chrom_df) == 0:
        max_chrom = 0

    if max_chrom > max_plas:
        return False
    elif max_plas > max_chrom:
        return True
    else:
        return False


def simplify_web(dict):
    for key in list(dict):
        dict[key].remove(key)
    simple=dict
    for key in list(dict):
        try:
            for element in dict[key]:
                if element in list(simple):
                    simple.pop(element)
        except:
            pass
    return simple

def return_contig(file_loc, contig_arr):
    for contig in contig_arr:
        if contig[0]==file_loc:
            return contig

def plasmid_writer(contig_info, tempdir):
    filepath=""
    for contig in contig_info:
        if contig[3]:
            filepath=f"{filepath} {str(contig[0])}"
    temp_plas_dir=f"{tempdir}/tot_plas_contigs"
    makedir(temp_plas_dir)
    temp_fas_path=f"{temp_plas_dir}/plasmid.fasta"
    cmd=f"cat {filepath} > {temp_fas_path}"
    proc1 = sp.run(cmd, shell=True, stdout=sp.PIPE)
    print(proc1.stdout)
    return temp_fas_path

def mash_triangle(plas_fas_path, plasdb, ext, tempdir, threads):
    concat_tri_dir=f"{tempdir}/triangle"

    makedir(concat_tri_dir)

    plasdb_path=f"{plasdb}{ext}"

    num_cont=len(read_fasta(plas_fas_path))

    temp_db_fas = f"{concat_tri_dir}/db.fasta"
    temp_plas_fas=f"{concat_tri_dir}/seq.fasta"

    temp_output=f"{concat_tri_dir}/dist.tsv"

    shutil.copy(plasdb_path, temp_db_fas)
    shutil.copy(plas_fas_path, temp_plas_fas)

    input_mash=f"{concat_tri_dir}/concat_tot.fasta"

    cmd_1=f"cat {temp_db_fas} {temp_plas_fas} > {input_mash}"

    cmd_2=f"mash sketch -i {input_mash} -p {threads} -o {input_mash}.msh"

    cmd_3=f"mash triangle -i {input_mash}.msh -p {threads} > {temp_output}"

    proc1 = sp.run(cmd_1, shell=True, stdout=sp.PIPE)

    proc2 = sp.run(cmd_2, shell=True, stdout=sp.PIPE)

    proc3 = sp.run(cmd_3, shell=True, stdout=sp.PIPE)

    print(cmd_1,cmd_2, cmd_3)

    print(proc1.stderr, proc2.stderr, proc3.stderr)

    return [num_cont, temp_output]

def cluster(dist_file, num_cont):
    with open(dist_file, 'r') as file:
        lines = file.readlines()
    lines=lines[1:]
    dimension_sq=len(lines[1:])
    arr = np.zeros((dimension_sq, dimension_sq), dtype=float)

    line_num = 0
    for line in lines:
        line=line.replace('\n','')
        line_arr=line.split("\t")[1:]
        entry_num=0

        for entry in line_arr:
            try:
                arr[line_num, entry_num]=entry
                try:
                    arr[entry_num, line_num]=entry
                except:
                    pass
                entry_num+=1
            except:
                pass
        line_num+=1
    print(arr)

    condensed_distance = squareform(arr)

    # Calculate the linkage matrix using single linkage
    linkage_matrix = linkage(condensed_distance, method='single')

    # Replace 'threshold_value' with your desired threshold for clustering
    threshold_value = 0.2

    # Get the clustered elements using fcluster
    clusters = fcluster(linkage_matrix, threshold_value, criterion='distance')
    results = clusters[-num_cont:]
    print("Clusters:", results)

    return results

def seq_arr_bplen(seq_arr):
    sum=0
    for sequence in seq_arr:
        sum+=len(sequence)
    return sum

def plasmid_splitter(results, plas_fas, location):
    plas_dict={}
    plas_fasta=read_fasta(plas_fas)
    for cluster in results:
        plas_dict[cluster]=[]

    for i in range(len(results)):
        cluster=results[i]
        plas=plas_fasta[i]

        plas_dict[cluster].append(plas)

    to_rem=[]

    for cluster in plas_dict:
        seq_arr=plas_dict[cluster]
        if seq_arr_bplen(seq_arr)<10000:
            to_rem.append(cluster)

    for rem in to_rem:
        plas_dict.pop(rem)

    for cluster in plas_dict:
        clus_name=f"{location}/plasmid_cluster_{cluster}.fasta"
        write_fasta(clus_name, plas_dict[cluster])

def main():
    args=arguments()
    tempdir=args.tempdir
    makedir(f"{tempdir}")

    contig_location = makedir(f"{tempdir}/contigs")
    blast_plas_location = makedir(f"{tempdir}/blast_plas_res")
    blast_chrom_location = makedir(f"{tempdir}/blast_chrom_res")

    fasta_breaker(args.infile, contig_location)
    contig_info=contig_class(contig_location, blast_plas_location, blast_chrom_location, args.plasdb, args.chromdb, args.num_threads)

    outputdir=args.outputdir
    makedir(outputdir)
    chrom_path=f"{outputdir}/chromosome.fasta"
    write_chrom_fasta(contig_info, chrom_path)
    plas_fas_path=plasmid_writer(contig_info, tempdir)
    mash_o=mash_triangle(plas_fas_path, args.plasdb, args.plasdb_ext,tempdir, args.num_threads)
    plas_results=cluster(mash_o[1], mash_o[0])
    plasmid_splitter(plas_results, plas_fas_path, outputdir)



    #contig_blaster(contig_location, blast_plas_location, blast_chrom_location, args.plasdb, args.chromdb, args.num_threads)


    shutil.rmtree(tempdir)



main()