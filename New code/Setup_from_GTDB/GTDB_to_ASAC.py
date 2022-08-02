"""
GTDB to ASAC
Returns a file res_accession.txt, containing the Assembly accession identifiers of all genomes within GTDB release

Usage:
    GTDB_to_ASAC <input_file> [-n <n>]
    

Options:
    -n=<n>      Name of resulting .txt file [default: res_accession.txt]
    -h --help   Show this screen.
"""

from docopt import docopt

def load_paths(ref_path):
    with open(ref_path.strip()) as f:
        lines = f.readlines()
        return(lines)

def GTDB_to_accession(GTDB_file):
    genomes = load_paths(GTDB_file)
    accession_dataset = list()
    for line in genomes:
        accession_dataset.append(line.split('\t')[0][3::])
    return(accession_dataset)

def print_ftp_file(ftp_filename, ftp_list):
    with open(ftp_filename, 'w') as f:
        for line in ftp_list:
            f.write(line+'\n')
    f.close()
    return

def main(args):
    GTDB_file = args['<input_file>']
    res_name = args['-n']
    accession_list = GTDB_to_accession(GTDB_file)
    print_ftp_file(res_name, accession_list)

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)