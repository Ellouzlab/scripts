import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
from utils import *

def arguments():
    parser = argparse.ArgumentParser(description="Finds Pantoea plasmids based on homology")

    parser.add_argument('-i', '--input', required=True, type=str, help="csv containing fasta")
    parser.add_argument('-f', '--fcol', required=True, type=str, help="Enter columns with genes e.g., 2,3,4 or 2-4")
    parser.add_argument('-n', '--ncol', required=True, type=str, help='column with what you want to be names e.g., 2,3,4 or 2-4')
    parser.add_argument('-o', '--output', required=False, default='output', type=str, help="output folder containing genes")
    parser.add_argument('-r', '--fill_missing', required=False, default='', type=str,help='Replace the length of missing genes with a character Default: nothing')
    args = parser.parse_args()
    return args

def is_empty_or_whitespace(s):
    return not s or s.isspace()


def get_col_arr(col_input):
    no_space_col=col_input.replace(' ','')
    parts = no_space_col.split(',')

    result = []

    for part in parts:
        if '-' in part:
            start, end = map(int, part.split('-'))
            result.extend(range(start, end + 1))
        else:
            result.append(int(part))


    return [x - 1 for x in result]

def add_n(df, f_arr, rep):
    for col in f_arr:
        col_seq=[]
        for row in df.index:
            value=str(df.iloc[row, col]).replace(' ','').split('\n')
            no_record_title=[]
            for line in value:
                if not '>' in line:
                    no_record_title.append(line)
            fasta_seq=''.join(no_record_title)
            df.iloc[row, col]=fasta_seq
            col_seq.append(fasta_seq)
        longest_length = max(len(s) for s in col_seq)

        for row in df.index:
            value=df.iloc[row,col]
            if is_empty_or_whitespace(value) or value=='nan':
                print(f"Row {row}, Column {col} replaced with {rep}")
                print(f"parsed as {value}")
                df.iloc[row, col]=rep*longest_length
            else:
                print(f"found sequence in Row {row}, Column {col}")
    return df

def df_to_dict(df, ncol, colnum):
    fasta_dict={}
    for row in df.index:
        id_arr=[]
        for col in ncol:
            id_arr.append(str(df.iloc[row, col]).replace(' ', '_'))
        fasta_id=('_'.join(id_arr)).replace(' ', '_')

        seq=str(df.iloc[row, colnum])
        seq=seq.replace('\n','')

        fasta_dict[fasta_id]=seq
    return fasta_dict


def main():
    args=arguments()
    input=args.input
    output=args.output
    ncol=args.ncol
    fcol=args.fcol
    rep=args.fill_missing

    n_arr=get_col_arr(ncol)
    f_arr=get_col_arr(fcol)

    makedir(output)

    tot_df=pd.read_csv(input)

    filled_df=add_n(tot_df, f_arr, rep)

    for col in f_arr:
        out_loc=f"{output}/column_{col}.fasta"
        sequences=df_to_dict(filled_df, n_arr, col)
        records = [SeqRecord(Seq(seq), id=seq_id, description='') for seq_id, seq in sequences.items()]
        write_fasta(out_loc, records)

main()
