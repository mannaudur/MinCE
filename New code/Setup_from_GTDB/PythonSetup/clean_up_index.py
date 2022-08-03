"""
Clean up Index
Accepts two index files (bacteria and archaea) and creates a large index file for MinCE

Usage:
    clean_up_index <bac_index> <arc_index> [-n <n>]
    

Options:
    -n=<n>      Name of resulting .txt file [default: res_accession.txt]
    -h --help   Show this screen.
"""

from docopt import docopt
import pandas as pd

def main(args):
    bac_index = pd.read_csv(args['<bac_index>'], sep="\t", names=["Id", "Accession", "GenBank_ftp", "RefSeq_ftp"])
    arc_index = pd.read_csv(args['<arc_index>'], sep="\t", names=["Id", "Accession", "GenBank_ftp", "RefSeq_ftp"])
    new_index = pd.concat([arc_index, bac_index], ignore_index=True)
    new_index.sort_values(by = "Id")
    reduced = new_index.loc[:,["Id","Accession"]]

    with open(args['-n'], 'w') as f:
        for i in range(len(reduced)):
            line = str(reduced.loc[i,"Id"])+"\t"+str(reduced.loc[i,"Accession"])+"\n"
            f.write(line)

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)