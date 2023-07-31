import os, sys, io, argparse
import pandas as pd
from codecs import decode
from Bio import SeqIO
import subprocess as sp

def main():
    '''Main function'''
    args = get_args()
    run(**args.__dict__)

def run(
        indir: str,
        outdir: str,
        perc_id: str,
        q_cov: str,
        e_val: str,
        query: str,
        mode: str,
        num_hsp:str

) -> None:
    tot_df=None
    txt=""
    res_txt=""
    i=0
    for prokka_dir in os.listdir(indir):
        prokka_dir_path=f"{indir}/{prokka_dir}"

        dir_path=db_dir_path(prokka_dir_path, mode)
        prokka=build_ProkkaDir(prokka_dir_path)

        if not os.path.exists(dir_path[0]):
            makedb_run(prokka, mode)

        output=blast_run(query, prokka, mode, perc_id, q_cov, num_hsp, e_val)
        result=parse_blast_result(output, mode)
        result=set_blast_df(prokka_dir, result)

        if i==0:
            tot_df=result
        else:
            tot_df=pd.concat([tot_df, result])
        i+=1

        if not "Empty DataFrame" in str(result):
            num_hits=len(result)
            txt+=f'{prokka_dir} had {num_hits} hits\n'
        else:
            txt+=f"{prokka_dir} had no hits\n"

    txt=txt+tot_df.to_string()
    print(txt)

    output=makedir(outdir, force_make=True)
    write_txt(f"{output}/output.txt", txt)
    print(output)


def build_ProkkaDir(path: str)-> object:
    '''
    build an object with file paths as locations
    :param path: path to prokka directory
    :return: an object that has a path to all important prokka objects
    '''

    class ProkkaDir:
        def __init__(self, ffn_file, fna_file, gff_file, log_file, faa_file, tbl_file,location):
            self.ffn = ffn_file
            self.fna = fna_file
            self.gff = gff_file
            self.log = log_file
            self.faa = faa_file
            self.tbl = tbl_file
            self.location=location

        def __call__(self):
            #return filename
            return self

    properties = {
        'ffn_file': None,
        'fna_file': None,
        'gff_file': None,
        'log_file': None,
        'faa_file': None,
        'tbl_file': None,
        'location': path
    }

    #get files by extension
    for file in os.listdir(path):
        file_path=path+"/"+file

        if file.split(".")[-1]=="ffn":
            properties['ffn_file']=file_path

        if file.split(".")[-1]=="fna":
            properties['fna_file']=file_path

        if file.split(".")[-1]=="gff":
            properties['gff_file']=file_path

        if file.split(".")[-1]=="log":
            properties['log_file']=file_path

        if file.split(".")[-1]=="faa":
            properties['faa_file']=file_path

        if file.split(".")[-1]=="tbl":
            properties['tbl_file']=file_path

    #print warning message if a filetype is not in prokkadir
    for ext in properties:
        if ext==None:
            print("WARNING: folder ", path, 'has no ', ext)

    #build ProkkaDir object
    ProkkaDir_obj=ProkkaDir(**properties)
    return ProkkaDir_obj

def set_blast_df(prokkadir, result):
    result["subject_seq_id"] = result["subject_title"].str.split(n=1).str[0]
    result["subject_title"] = result["subject_title"].str.split(n=1).str[1]
    result["sample"] = prokkadir
    result = set_column_sequence(result, ["sample"])
    return result

def blast_run(query, prokka, mode: str, perc, qcov, num_hsp, eval):
    db_path=db_dir_path(prokka.location, mode)[0]
    name=db_dir_path(prokka.location, mode)[1]
    if mode=='nucl':
        cmd=f"blastn -query {query} -perc_identity {perc} -evalue {eval} -max_hsps {num_hsp} -qcov_hsp_perc {qcov} " \
            f"-outfmt '10 stitle qlen slen length qcovhsp pident evalue' -db {db_path}/{name}"
    else:
        cmd=f"blastn -query {query} -perc_identity {perc} -evalue {eval} -max_hsps {num_hsp} -qcov_hsp_perc {qcov} " \
            f"-outfmt '10 stitle qlen slen length qcovhsp pident evalue' -db {db_path}/{name}"
    proc1=sp.run(cmd, shell=True, stdout=sp.PIPE)
    return proc1.stdout

def parse_blast_result(result, mode):
    if mode=='nucl':
        col_names=["subject_title", "query_length", "subject_length", "alignment_length", "query_coverage", "percent_identical", "e value"]
        df=pd.read_csv(
            io.StringIO(decode(result)),
            names=col_names
        )
        return df

def makedb_run(prokka, mode: str):
    db_path=db_dir_path(prokka.location, mode)[0]
    name=db_dir_path(prokka.location, mode)[1]
    makedir(db_path)
    print(db_path)
    if mode=='nucl':
        cmd=f"makeblastdb -in {prokka.ffn} -dbtype {mode} -out {db_path}/{name}"
    else:
        cmd = f"makeblastdb -in {prokka.faa} -dbtype {mode} -out {db_path}/{name}"
    sp.run(cmd)


def db_dir_path(seq_path: str, mode: str)->list:
    basename = seq_path.split('/')[-1]
    no_ext = basename.split('.')[0]

    dir_back = seq_path+"/../.."
    db_dir=f"{dir_back}/db_{mode}/{no_ext}"

    return [db_dir, no_ext]

def makedir(path: str, force_make=False) -> str:
    '''
    Makes directory for the given path
    :param path: the location of results
    :param force_make: forces a directory to be made for results with a slight change in folder name
    :return: returns path to new directory
    '''

    #try making a directory
    try:
        os.mkdirs(path)
        print('worked')
        return path
    except:
        print('failed')
        if not force_make:
            return path

    #forces a directory to be made by adding a number at the end
    i = 1
    if force_make:
        while True:
            new_name = path + '(' + str(i)+')'
            try:
                os.mkdirs(new_name)
                break
            except:
                i+=1
        return new_name

def write_fasta(path: str, recordlist: list)->None:
    '''
    writes a fasta file
    :param location of fastafile:
    :param sequences in fastafile:
    :return: None
    '''
    with open(path, "w") as output_handle:
        SeqIO.write(recordlist, output_handle, "fasta")

def write_txt(path: str, text:str)->None:
    '''
    writes a txt
    :param path: location of file
    :param dictionary: results
    :return: None
    '''
    with open(path, "w") as f:
        f.write(text)

def set_column_sequence(dataframe, seq, front=True):
    '''
    Takes a dataframe and a subsequence of its columns,
    returns dataframe with seq as first columns if "front" is True,
    and seq as last columns if "front" is False.
    '''
    cols = seq[:] # copy so we don't mutate seq
    for x in dataframe.columns:
        if x not in cols:
            if front: #we want "seq" to be in the front
                #so append current column to the end of the list
                cols.append(x)
            else:
                #we want "seq" to be last, so insert this
                #column in the front of the new column list
                #"cols" we are building:
                cols.insert(0, x)
    return dataframe[cols]

def get_args() -> argparse.Namespace:
    """
    gets arguments
    :return: arguments
    """
    description = "Searches prokka output directories and returns a fastafile and summary file"
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        type=str,
        help="Directory containing folders with prokka output directories",
        metavar="INDIR"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Output directory",
        metavar="OUTDIR",
        default="BlastFind_Results"
    )
    parser.add_argument(
        "-p",
        "--perc_id",
        type=str,
        help="min percent identity (Default: 50)",
        metavar="perc_id",
        default='50'
    )
    parser.add_argument(
        "-c",
        "--q_cov",
        type=str,
        help="min query cover (Default: 20)",
        metavar="perc_id",
        default='20'
    )
    parser.add_argument(
        "-e",
        "--e_val",
        type=str,
        help="min e-value",
        metavar="e-value",
        default="1E06"
    )
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        required=True,
        help="Search for this fasta",
        metavar="query_string",
    )
    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        type=str,
        help="Mode of blasting: 'nucl' (DNA) or 'prot' (protein) ",
        choices=['nucl', 'prot'],
        metavar="nucl/prot",
        default="prot"
    )
    parser.add_argument(
        "-n",
        "--num_hsp",
        type=str,
        help="maximum number of hits per genome file (Default: 5)",
        metavar="max_hsps",
        default=5
    )
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()
