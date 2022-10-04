"""
Translate

Usage:
    translate <input_file>

Options:
    -h --help   Show this screen.
"""

import json
import gzip
from docopt import docopt



def load_paths(ref_path):
    with open(ref_path.strip()) as f:
        lines = f.readlines()
        lines = [x.strip() for x in lines]
        return(lines)
    
def parse_fasta(fn):
    desc = []
    seqs = []

    s = []

    with gzip.open(fn,mode='rt') as fin:        
        for line in fin:
            if line[0] == '>':
                desc.append(line.strip()[1:])
                if s:
                    seqs.append(''.join(s))
                s = []
            else:
                s.append(line.strip())

        seqs.append(''.join(s))

    return {'file': fn, 'desc': desc, 'seqs': seqs}

def main(args):
    input_file = args['<input_file>']

    pseudonyms = dict()

    for path in load_paths(input_file):
        subnames = [x.split(' ')[0] for x in parse_fasta(path)['desc']]
        proper_name = path.split('/')[-1].split('_genomic.fna.gz')[0]
        for name in subnames:
            pseudonyms[name] = proper_name

    with open(input_file+".json", "w") as outfile:
        json.dump(pseudonyms, outfile)


if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)