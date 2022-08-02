"""
Retry missing genomes
Returns a file res_rsync.txt. Start downloading files by running:
rsync rsync --copy-links -va --no-relative rsync://ftp.ncbi.nlm.nih.gov/genomes/all/ --files-from=res_rsync.txt archaea_207/

Usage:
    retry_missing_genomes <found_file> <index_file> [-n <n>]
    

Options:
    -n=<n>      Name of resulting .txt file [default: retry_missing.txt]
    -h --help   Show this screen.
"""

from docopt import docopt
import pandas as pd
import re

def load_paths(ref_path):
    with open(ref_path.strip()) as f:
        lines = f.readlines()
        return(lines)

def acc_to_path(csv_line):
    acc = csv_line[1]
    ass = str(csv_line[2]).replace(" ", "_")
    fas = acc+'_'+ass
    match = re.search('(\w+)_(\d+)\.(\d)', acc)
    if match:
        prefix = match.group(1)
        accn = match.group(2)
        accp = re.findall('.{3}', accn)
        accp = '/'.join(accp)
        return (prefix, accp, fas, fas+'_genomic.fna.gz')
    else:
        raise ValueError("could not get FTP path from ".format(acc))

def main(args):
    found_file = args['<found_file>']
    index_file = args['<index_file>']
    outfile = args['-n']
    index_table = pd.read_csv(index_file, sep="\t", header=None)
    found = load_paths(found_file)
    found_set = set()
    for mem in found:
        found_set.add(mem.split('/')[-1][0:15])

    wanting = index_table[index_table.iloc[:,1].isin(found_set) == False]

    for i in range(len(wanting[1])):
        if wanting.iloc[i,1][2] == 'F':
            wanting.iloc[i,1] = 'GCA'+wanting.iloc[i,1][3:]
        elif wanting.iloc[i,1][2] == 'A':
            wanting.iloc[i,1] = 'GCF'+wanting.iloc[i,1][3:]
        else:
            print('Member '+mem+' is found to be weird')

    if(len(wanting) == 0):
        print('No genomes missing')

    with open(outfile, 'w') as f:
        for i in range(len(wanting)):
            f.write(('/').join(acc_to_path(wanting.iloc[i]))+"\n")

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)