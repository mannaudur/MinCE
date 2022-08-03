"""
Index to RSYNC
Returns a file res_rsync.txt. Start downloading files by running:
rsync rsync --copy-links -va --no-relative rsync://ftp.ncbi.nlm.nih.gov/genomes/all/ --files-from=res_rsync.txt archaea_207/

Usage:
    Index_to_RSYNC <index_file> [-r <r>] [-m <m>]
    

Options:
    -r=<r>      Name of resulting .txt file [default: res_rsync.txt]
    -m=<m>      Name of file for missing genomes [default: res_missing.txt]
    -h --help   Show this screen.
"""

from docopt import docopt
import pandas as pd



def main(args):
    index_file = args['<index_file>']
    rsync_found = args['-r']
    rsync_missing = args['-m']

    ftp_list = list()
    missing = list()

    index_file = pd.read_csv(index_file, sep="\t", names=["Id", "Accession", "GenBank_ftp", "RefSeq_ftp"])
    for i in range(len(index_file)):
        link = index_file.loc[i,"GenBank_ftp"]
        if pd.isna(link):
            link = index_file.loc[i,"RefSeq_ftp"]
            if pd.isna(link):
                name = index_file.loc[i,"Accession"].replace('GCA', '%temp%').replace('GCF', 'GCA').replace('%temp%', 'GCF')
                missing.append(name+'\n')
                continue
        prefix = link.split('/genomes/all/')[-1]
        name = prefix.split('/')[-1]+'_genomic.fna.gz\n'
        ftp_list.append(prefix+'/'+name)

    with open(rsync_found, 'w') as f:
        for line in ftp_list:
            f.write(line)

    with open(rsync_missing, 'w') as f:
        for line in missing:
            f.write(line)
    

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)