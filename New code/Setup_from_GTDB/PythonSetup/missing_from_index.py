"""
Missing from Index
Accepts two an index file and an accession file, returns files present in accession file but missing from index file.

Usage:
    missing_from_index <index_file> <accession_file> [-n <n>]
    

Options:
    -n=<n>      Name of resulting .txt file [default: accession_missing.txt]
    -h --help   Show this screen.
"""

from docopt import docopt
import pandas as pd

def strip_list(lines):
    for i in range(len(lines)):
        lines[i] = lines[i].strip()
    return(lines)

def load_paths(ref_path):
    with open(ref_path.strip()) as f:
        lines = f.readlines()
        return(lines)

def main(args):
    viral = pd.read_csv(args['<index_file>'], sep="\t", names=["Id", "Accession", "GenBank_ftp", "RefSeq_ftp"])
    org = strip_list(load_paths(args['<accession_file>']))

    found = set(viral.loc[:,"Accession"])

    missing = set(org) - found
    with open(args['-n'], 'w') as f:
        for mem in missing:
            f.write(mem+'\n')


if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)