"""
Index to RSYNC
Returns a file res_rsync.txt. Start downloading files by running:
rsync rsync --copy-links -va --no-relative rsync://ftp.ncbi.nlm.nih.gov/genomes/all/ --files-from=res_rsync.txt archaea_207/

Usage:
    Index_to_RSYNC <input_file> [-n <n>]
    

Options:
    -n=<n>      Name of resulting .txt file [default: res_rsync.txt]
    -h --help   Show this screen.
"""

from docopt import docopt
import pandas as pd
import re

# Taken (and amended) from https://github.com/ctSkennerton/scriptShed/blob/master/download_ncbi_assembly.py
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
    index_file = args['<input_file>']
    outfile = args['-n']
    ar53 = pd.read_csv(index_file, sep="\t", header=None)
    with open(outfile, 'w') as f:
        for i in range(len(ar53)):
            f.write(('/').join(acc_to_path(ar53.iloc[i]))+"\n")

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)