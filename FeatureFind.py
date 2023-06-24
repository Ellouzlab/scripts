import os
from Bio import SeqIO
import argparse

def main():
    '''Main function'''
    args = get_args()
    run(**args.__dict__)

def run(
        indir: str,
        outdir: str,
        caps: str,
        excl: str,
        query: str,
        output_type: str
) -> None:
    '''
    run this program
    '''

    #handles incorrect input
    if not (caps=='y' or caps=='n'):
        exit(print("-c/--caps only accepts the option 'y' or 'n'"))

    if not (output_type=="DNA" or output_type=="Protein"):
        exit(print("-t/--output_type only accepts the options 'Protein' or 'DNA'"+ output_type))

    #get absolute paths
    outdir_true=os.getcwd()+"/"+outdir
    indir_true=os.getcwd()+"/"+indir

    #build directory for outputs and save their location
    dir=makedir(outdir_true)
    print(dir)
    query_outputdir=outdir_true+"/"+query
    path_to_dir=makedir(query_outputdir, force_make=True)


    #build a ProkkaDir object for each prokka output directory and appends to array
    ProkkaDirArray=[]
    for folder in os.listdir(indir_true):
        folder_path=indir_true+"/"+folder
        ProkkaDirArray.append(
            build_ProkkaDir(folder_path)
        )

    final_res=dict()
    combined_fasta=[]
    for ProkkaDir in ProkkaDirArray:
        if output_type=="DNA":
            recordlist=file_reader(ProkkaDir.ffn)
        else:
            recordlist=file_reader(ProkkaDir.faa)

        results=gene_find(query, excl, recordlist, caps)

        for sequence in results:
            combined_fasta.append(sequence)

        final_res[ProkkaDir()]=len(results)

        print(results)
        print(final_res)
    write_txt(path_to_dir+"/"+"Query_result.txt", final_res)
    write_fasta(path_to_dir+"/"+"Query_result.fasta", combined_fasta)

def get_args() -> argparse.Namespace:
    """
    gets arguments
    :return: arguments
    """
    description = "Searches prokka output directories and returns a fasta and summary file"
    parser = argparse.ArgumentParser(description=description, add_help=True)
    parser.add_argument(
        "-i",
        required=True,
        type=str,
        help="Directory containing folders with prokka output directories",
        metavar="INDIR"
    )
    parser.add_argument(
        "-o",
        type=str,
        help="Output directory",
        metavar="OUTDIR",
        default="FeatureFind_Results"
    )
    parser.add_argument(
        "-c",
        type=str,
        help="ignore caps (y/n) (Default: y)",
        metavar="'y'/'n'",
        default='y'
    )
    parser.add_argument(
        "-e",
        type=str,
        help="No suggestions will be provided if their description contain this string",
        metavar="exclusion_string"
    )
    parser.add_argument(
        "-q",
        required=True,
        type=str,
        help="Search for this string in record descriptions",
        metavar="query_string"
    )
    parser.add_argument(
        "-t",
        type=str,
        help="type of outputted fasta (Protein or DNA) (Default: DNA)",
        metavar="'Protein'/'DNA'",
        default="DNA"
    )
    args = parser.parse_args()
    return args

def build_ProkkaDir(path: str)-> object:
    '''
    build an object with file paths as locations
    :param path: path to prokka directory
    :return: an object that has a path to all important prokka objects
    '''

    class ProkkaDir:
        def __init__(self, ffn_file, fna_file, gff_file, log_file, faa_file, tbl_file):
            self.ffn = ffn_file
            self.fna = fna_file
            self.gff = gff_file
            self.log = log_file
            self.faa = faa_file
            self.tbl = tbl_file

        def __call__(self):
            #return filename
            return self.ffn.split('/')[-1].split('.')[0]

    properties = {
        'ffn_file': None,
        'fna_file': None,
        'gff_file': None,
        'log_file': None,
        'faa_file': None,
        'tbl_file': None
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


def makedir(path: str, force_make=False) -> str:
    '''
    Makes directory for the given path
    :param path: the location of results
    :param force_make: forces a directory to be made for results with a slight change in folder name
    :return: returns path to new directory
    '''

    #try making a directory
    try:
        os.mkdir(path)
        return path
    except:
        if not force_make:
            return path

    #forces a directory to be made by adding a number at the end
    i = 1
    if force_make:
        while True:
            new_name = path + '(' + str(i)+')'
            try:
                os.mkdir(new_name)
                break
            except:
                i+=1
        return new_name



def file_reader(fastafile: str, changeid=True) -> list:
    '''
    Reads file and returns sequence objects in an array
    :param fastafile: path to file
    :return: list of sequence objects
    '''
    recordlist = []
    #open file
    with open(fastafile) as handle:
        try:
            for record in SeqIO.parse(handle, "fasta"):
                if changeid:
                    record.id=fastafile.split('/')[-1].split('.')[0]
                recordlist.append(record)
        except:
            print('could not read the sequence file ', fastafile)
    return recordlist

def write_fasta(path: str, recordlist: list)->None:
    '''
    writes a fasta file
    :param location of fastafile:
    :param sequences in fastafile:
    :return: None
    '''
    with open(path, "w") as output_handle:
        SeqIO.write(recordlist, output_handle, "fasta")

def write_txt(path: str, dictionary: dict)->None:
    '''
    writes a txt
    :param path: location of file
    :param dictionary: results
    :return: None
    '''

    message=""
    for element in dictionary:
        message+=str(element)+" has "+str(dictionary[element])+" hits\n"
    print(message)

    with open(path, "w") as f:
        f.write(message)
def gene_find(query: str, exclusion: str, recordlist: list, caps:str)->list:
    '''
    Returns records that match criteria
    :param query: The string being searched for in the descriptions
    :param exclusion: If this string is contained within a result, it is excluded
    :param recordlist: Array of sequence objects to be searched
    :param caps: to ignore or not to ignore capitalizations
    :return: list of sequence objects that match criteria
    '''

    #if no exclusion, make the exclusion string something riduculous so nothing is excluded
    if exclusion==None:
        #a random sequence of charachters that won't be found in descriptions
        exclusion="AUUIEUIUEIYWIEUIUIEUIUWIUwMNJCJHDJF#%@!%!!#$%"

    #Search record descriptions
    hits=[]
    for record in recordlist:
        if caps=='y':
            if query.lower() in record.description.lower() and not exclusion.lower() in record.description.lower():
                hits.append(record)
        else:
            if query in record.description and not exclusion in record.description:
                hits.append(record)
    return hits


if __name__ == "__main__":
    main()