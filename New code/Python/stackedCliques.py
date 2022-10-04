"""
StackedCliques

Usage:
    stackedCliques <input_file> [-t <t>]

Options:
    -t=<t>      Largest sketch distance value to find distinguishing sequences [default: 3]
    -h --help   Show this screen.
"""

from re import I
from docopt import docopt
import numpy as np
import gc

def load_sketch(path):
    fixed_path = 'sketches/'+path.split('/')[-1]+'.sketch'
    with open(fixed_path.strip()) as f:
        lines = f.readlines()
        name = lines[0]
        sketch = set(lines[4::])
        return(sketch, name)
    
def load_paths(ref_path):
    with open(ref_path.strip()) as f:
        lines = f.readlines()
        return(lines)
    
def distMat(cluster_path):
    paths = load_paths(cluster_path)
    n = len(paths)
    D = np.zeros((n,n))
    names = []
    for i in range(n):
        sketch, name = load_sketch(paths[i].strip())
        names.append(name.strip())
        for j in range(n):
            D[i,j] = 5000 - len(sketch.intersection(load_sketch(paths[j].strip())[0]))
    return(D, names)
    
def simple_cliques_with_subset_merge(D, thresh):
    cliques = list()
    mems = set(np.where(D == thresh)[0])
    for mem in mems: #members in the matrix
        merge = set(np.where(D[mem,:] == thresh)[0]) #Using ids to get original indices for accounting
        if(len(merge)):
            add_ = True
            for clq in cliques:
                if clq.issubset(merge) or merge.issubset(clq): #Check if we have a superset in the log
                    clq.update(merge) #merge subset with superset
                    add_ = False #don't make a new one
            if add_:
                cliques.append(merge) #add a new set

    for i in range(len(cliques)): #Second pass to join cliques differing by 1 member
        for j in range(i,len(cliques)): #Iterating over "downstream" cliques
            smaller_clique = min(len(cliques[i]), len(cliques[j])) #Size of smaller clique
            if len(cliques[i].intersection(cliques[j])) - 1 <= smaller_clique:
                cliques[i].update(cliques[j])
                cliques[j] = set()

    cliques = [x for x in cliques if len(x)>0]

    return(cliques)
    
def stacked_cliques(D, thresh):
    cliques = []
        
    for t in range(thresh+1):
        clq = simple_cliques_with_subset_merge(D,t)
        if(len(clq) > 0):
            cliques.append([t, clq])
    return(cliques)

def write_cliques_to_disk(clique_log, clique_list, names, cluster_path):
    cluster_id = cluster_path.split('/')[-1].split('.')[0]
    inner_cluster_dict = dict()
    for i in range(len(clique_list)):
        thresh = clique_list[i][0]
        cliques = clique_list[i][1]
        for j in range(len(cliques)):
            members = cliques[j]
            clique_id = cluster_id+"_"+str(thresh)+"_"+str(j)
            clique_log.append([clique_id, len(members)])
            with open('cliques/'+clique_id+".clique", 'w') as f:
                for mem in members:
                    f.write(names[mem]+'\n')
                    if mem in inner_cluster_dict:
                        inner_cluster_dict[mem].append(clique_id)
                    else:
                        inner_cluster_dict[mem] = [clique_id]
    return(clique_log)

def print_clique_log(clique_log, thresh):
    clique_log.sort(key = lambda x: x[1], reverse = True)
    with open('Overview_of_cliques', 'w') as f:
        f.write("Cliques formed by iterating up to and including sketch distance "+str(thresh)+"\n")
        f.write("Number of cliques: ")
        f.write(str(len(clique_log)) + "\n\n")
        for line in clique_log:
            f.write((" ").join(str(entry) for entry in line))
            f.write('\n')

def print_new_index(indices):
    with open('index', 'w') as f:
        for line in indices:
            f.write(line.strip())
            f.write('\n')
                        
def clique_every_cluster(all_file, thresh):
    all_paths = load_paths(all_file)
    clique_log = []
    i = 0
    for cluster_path in all_paths:
        D, names = distMat(cluster_path)
        clique_list = stacked_cliques(D, thresh)
        clique_log = write_cliques_to_disk(clique_log, clique_list, names, cluster_path)
        i += 1
        print(i, "/", len(all_paths), "done")
        gc.collect
    print_clique_log(clique_log, thresh)
    
def main(args):
    all_file = args['<input_file>']
    thresh = int(args['-t'])
    clique_every_cluster(all_file, thresh)
    print("All done!")


if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)