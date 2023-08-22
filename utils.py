from Bio import SeqIO
import os, io
import pandas as pd
from codecs import decode
import subprocess as sp

def read_fasta(fastafile: str)->list:
    '''
    Reads fastafile
    :param fastafile: file path
    :return: list of fasta record
    '''
    recordlist = []
    # open file
    with open(fastafile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            recordlist.append(record)
    return recordlist


def write_fasta(path: str, recordlist: list)->None:
    '''
    Writes a fasta file to a given location
    :param path: location to write fasta
    :param recordlist: list of fasta records
    :return: None
    '''
    with open(path, "w") as output_handle:
        SeqIO.write(recordlist, output_handle, "fasta")


def makedir(path: str, force_make: bool=False)->str:
    '''
    Makes a directory if the given direcotry doesn't exist. If force_make true, makes a directory with a number
    :param path: location of new directory
    :param force_make: To make a new directory with a number if a directory already exists at given path
    :return: path to new dir
    '''
    try:
        os.mkdir(path)
        return path
    except:
        if not force_make:
            return path
    i = 1
    if force_make:
        while True:
            new_name = f"{path}({i})"
            try:
                os.mkdir(new_name)
                break
            except:
                i += 1
        return new_name


def fasta_breaker(fasta: list, tempdir: str)-> None:
    '''
    Writes a fasta list as individual files under contig 1 to n number of records.
    :param fasta: fasta records list
    :param tempdir: directory to write records to
    :return: None
    '''
    records = read_fasta(fasta)
    i = 1
    for contig in records:
        fasta_location = f"{tempdir}/contig_{i}"
        write_fasta(fasta_location, contig)
        i += 1

def read_fasta_single(fastafile: str)->object:
    '''
    Returns a single fasta record (only use for fasta with a single record)
    :param fastafile: path to fastafile
    :return: return fasta record
    '''
    with open(fastafile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return record

def blast2df(result: object)->object:
    '''
    Takes output from blast_run and converts to a dataframe
    :param result: output from blast_run
    :return: dataframe
    '''
    col_names=["qtitle", "query_length", "subject_length", "alignment_length", "query_coverage", "percent_identical", "e value", "qstart", "qend"]
    df=pd.read_csv(
        io.StringIO(decode(result)),
        names=col_names
    )
    df["Match_base_pid"]=df["query_coverage"]*df["percent_identical"]/100
    df=df.drop_duplicates(subset="qtitle")
    return df

def blast_run(query: str, db_path: str, perc:float=80.0, eval: float=0.00001, threads: int=1)->object:
    '''
    Runs blastn and returns a dataframe
    :param query: path to query fasta
    :param db_path: path to database for blasting
    :param perc: minimum percent id
    :param eval: minumum e-value
    :param threads: number of threads to use
    :return: dataframe
    '''
    cmd=f"blastn -query {query} -perc_identity {perc} -evalue {eval} -num_threads {threads} " \
        f"-outfmt '10 sseqid qlen slen length qcovus pident evalue qstart qend' -db {db_path}"
    proc1=sp.run(cmd, shell=True, stdout=sp.PIPE)

    df=blast2df(proc1.stdout)
    return df



def write_txt(path: str, text: str)-> None:
    '''
    writes a string
    :param path: path to txt
    :param text: string to write
    :return: None
    '''
    with open(path, "w") as f:
        f.write(text)

def plasclust2fasta(contig_arr, outdir):
    rep_arr=[]
    for contig in contig_arr:
        if contig[4]:
            rep_arr.append(contig[6])
    name='plasmid_'
    for rep in rep_arr:
        name+=str(rep)+'_'

    seq_arr=[]
    for contig in contig_arr:
        seq_arr.append(read_fasta_single(contig[0]))

    location=f"{outdir}/{name}.fasta"
    write_fasta(location, seq_arr)

def unclust2fasta(contig_arr, outdir):
    filepath=f"{outdir}/unclustered_plas_contigs.fasta"
    seq_arr=[]
    for contig in contig_arr:
        seq_arr.append(read_fasta_single(contig[0]))

    write_fasta(filepath, seq_arr)

