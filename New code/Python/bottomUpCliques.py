"""
BottomUpCliques

Usage:
    bottomUpCliques <input_file> [-t <t>] [-i <i>]

Options:
    -t=<t>      How big an atom has to be to be broken up (bigger than t) [default: 5]
    -i=<i>      At which distance to start joining singletons to nearest cliques [default: 3]
    -h --help   Show this screen.
"""

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
    
def read_indices_to_rev_dict(path):
    indices = load_paths(path)
    rev_indices = dict()
    for line in indices:
        entries = line.split(" ")
        rev_indices[entries[1].strip()] = entries[0]
    return(rev_indices)
    
def distMat(cluster_path, rev_indices):
    paths = load_paths(cluster_path)
    n = len(paths)
    D = np.zeros((n,n))
    names = []
    ids = []
    for i in range(n):
        sketch, name = load_sketch(paths[i].strip())
        names.append(name.strip())
        ids.append(rev_indices[name.strip()])
        for j in range(n):
            D[i,j] = 5000 - len(sketch.intersection(load_sketch(paths[j].strip())[0]))
    return(D, names, ids)

# Takes in a distance matrix D_, which starts out as the original distance matrix
# but with each iteration we work with a submatrix of D, corresponding to the
# sketches yet to be assigned to a sub-clique. The parameter ids is a mapping between
# the indices of the orignal distance matrix D and the new submatrix D_.
# Thresh refers to the threshold distance value to unite sketches into cliques
# and the reference matrix is a table showing where each sketch ended up, which
# is persistently appended throughout the iterative process (that is between runs
# of this algorithm)
def simple_cliques_with_subset_merge(D, ids, thresh, reference_matrix):
    mems = np.arange(D.shape[0]) #members in the submatrix
    cliques = list()
    singles = set()
    for mem in mems:
        merge = ids[np.where(D[mem,:] <= thresh)[0]] #Using ids to get original indices for accounting
        if len(merge)>1: #checking if there is any match within the designated threshold
            merge = set(merge)
            add_ = True
            for clq in cliques:
                #Check if we have a superset in the log or whether difference is only 1 sketch
                if clq.issubset(merge) or merge.issubset(clq) or min(len(clq-merge),len(merge-clq) < 2):
                    clq.update(merge) #merge subset with superset
                    add_ = False #don't make a new one
            if add_:
                cliques.append(merge) #add a new set
        else:
            if(len(merge)):
                singles.add(merge[0]) #adding to singletons if no match
    #print("For threshold", thresh, ":")
    #print(len(singles), " singletons found.")
    #print(len(cliques), " cliques found")
    for i, clique in enumerate(cliques):
        #print("Size:", len(clique), "-", clique)
        for mem in clique:
            #log the position of each member for reference
            #thresh refers to line# in ref_matrix and i refers to # in list of cliques for that threshold
            #for example, if ref_mat[k,j] = m, then sketch j was joined to clique m with threshold k.
            reference_matrix[thresh, mem] = i 
    return(cliques, singles, reference_matrix)

# Organizing step in between runs of simple_cliques_with_subset_merge()
def clean_up_cliques(D, ids_):
    ids = np.array(list(ids_)) #convert singletons to new array of ids
    D_ = D[np.ix_(ids,ids)] #create submatrix of original distance matrix D, pertaining to those singleton ids
    return(D_, ids)

# Runs iterative calls on simple_cliques_with_subset_merge() with ever larger threshold values
# up to and including iter_. The path is just the path to give to distMat. The clean_up paramter
# is a boolean designating whether to return singletons that couldn't be placed in a clique in 
# the given amount of iter_ (clean_up=False) or to iterate over left over singletons and the end
# of iteration and assign them to the cliques containing their closest pairing in the distance matrix
def bottom_up_cliques(D, iter_, clean_up = True):
    reference_matrix = np.ones((iter_+2,D.shape[0]))*-1 #Reference matrix with extra line for unassigned singletons
    cliques = []
    #clq is threefold: clq[0] is cliques found, clq[1] is singletons and clq[2] is updated reference matrix
    clq = simple_cliques_with_subset_merge(D,np.arange(D.shape[0]),0, reference_matrix)
    # cliques.append(clq[0:2]) #If we wanted singletons as well
    cliques.append(clq[0]) #Stowing cliques away for reference
    reference_matrix = clq[2] #Filling in the reference matrix

    if len(clq[1]) == 0: #if no returned singletons, the algorithm has finished
        #print("--- Reached completion at threshold iteration",i,"---")
        return(cliques)
        
    for i in range(iter_):
        D_,ids = clean_up_cliques(D, clq[1]) #Clean up stage between runs
        clq = simple_cliques_with_subset_merge(D_,ids,i+1, reference_matrix)
        # cliques.append(clq[0:2]) #If we wanted singletons as well
        cliques.append(clq[0])
        reference_matrix = clq[2]
        if len(clq[1]) == 0: #if no returned singletons, the algorithm has finished
            #print("--- Reached completion at threshold iteration",i,"---")
            clean_up = False
            break
    singletons = ids #Just a renaming for clarifying purposes
    if(clean_up):
        for singleton in singletons: #for each singleton left
            D[singleton,singleton] = float('inf') #we wanna find the actual closest match
            closest_dist = np.min(D[singleton,:])
            closest_sketches = np.where(D[singleton,:] == closest_dist)[0] #which tells us which member to look for
            for close_sketch in closest_sketches: #We're gonna add this to every clique in mininum distance - effectively overlap
                clique_thresh = np.argmax(reference_matrix[:,close_sketch]) #which was put into a clique with this threshold
                clique_id = reference_matrix[clique_thresh,close_sketch].astype(np.int64) #and this is its location in the list
                if clique_id < 0: # there's a chance the close_sketch is a singleton, though it's somewhat unlikely
                    cliques.append([set([singleton, close_sketch])])
                    clique_thresh = min(D[singleton, close_sketch].astype(int), D.shape[0]-1)
                    reference_matrix[clique_thresh,singleton] = clique_id
                    reference_matrix[clique_thresh,close_sketch] = clique_id
                else:
                    cliques[clique_thresh][clique_id].add(singleton) #so we add that singleton to its closest neighbour's clique
                    reference_matrix[clique_thresh,singleton] = clique_id #and log it so placed
                #print("Added singleton", singleton,"to clique no.",clique_id,"with threshold",clique_thresh)
                #print("Amended form of clique is", cliques[clique_thresh][clique_id])   
    else:
        for singleton in singletons:
            reference_matrix[iter_+1,singleton] = -2 #if clean_up = False, just log those singletons as unassigned
                
            
    return(cliques)

def write_cliques_to_disk(new_index, clique_log, cliques_, names, ids, cluster_path):
    cluster_id = cluster_path.split('/')[-1].split('.')[0]
    inner_cluster_dict = dict()
    for i in range(len(cliques_)):
        if(cliques_[i]):
            for j in range(len(cliques_[i])):
                members = cliques_[i][j]
                clique_id = cluster_id+"_"+str(i)+str(j)
                clique_log.append([clique_id, len(members)])
                with open('cliques/'+clique_id+".clique", 'w') as f:
                    for mem in members:
                        f.write(names[mem]+'\n')
                        if mem in inner_cluster_dict:
                            inner_cluster_dict[mem].append(clique_id)
                        else:
                            inner_cluster_dict[mem] = [clique_id]
    for mem in inner_cluster_dict:
        clique_id_expanded = (";").join(inner_cluster_dict[mem])
        new_index.append([ids[mem], names[mem], cluster_id, clique_id_expanded])
    return(new_index, clique_log)

def write_unamended_clique_to_disk(new_index, clique_log, names, ids, cluster_path):
    cluster_id = cluster_path.split('/')[-1].split('.')[0]
    clique_id = cluster_id+"_0"
    clique_log.append([clique_id, len(names)])
    with open('cliques/'+clique_id+".clique", 'w') as f:
        for i,name in enumerate(names):
            f.write(name+'\n')
            new_index.append([ids[i], name, clique_id])
    return(new_index, clique_log)

def print_clique_log(clique_log, iter_, thresh):
    clique_log.sort(key = lambda x: x[1], reverse = True)
    with open('Overview_of_cliques', 'w') as f:
        f.write("Cliques formed with parameters: \n")
        f.write("Breaking up atoms bigger than "+str(thresh)+"\n")
        f.write("Iterating bottoms up until distance "+str(iter_)+"\n")
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
                        
def clique_every_cluster(all_file, iter_, thresh):
    all_paths = load_paths(all_file)
    rev_indices = read_indices_to_rev_dict('indices')
    new_index = []
    clique_log = []
    i = 0
    for cluster_path in all_paths:
        D, names, ids = distMat(cluster_path, rev_indices)
        if len(names) > thresh:
            cliques_ = bottom_up_cliques(D, iter_, True)
            new_index, clique_log = write_cliques_to_disk(new_index, clique_log, cliques_, names, ids, cluster_path)
        else:
            new_index, clique_log = write_unamended_clique_to_disk(new_index, clique_log, names, ids, cluster_path)
        i += 1
        print(i, "/", len(all_paths), "done")
        gc.collect
    del rev_indices
    gc.collect
    print_clique_log(clique_log, iter_, thresh)
    indices = load_paths('indices')
    for line in new_index:
        indices[int(line[0])] = (" ").join(str(entry) for entry in line)+'\n'
    print_new_index(indices)
    
def main(args):
    all_file = args['<input_file>']
    threshold = int(args['-t'])
    iter_ = int(args['-i'])-1
    clique_every_cluster(all_file, iter_, threshold)
    print("All done!")


if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)