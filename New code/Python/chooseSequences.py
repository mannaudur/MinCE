"""
ChooseSequences

Usage:
    chooseSequences <input_file> [-k <k>] [-n <n>] [-B <B>]

Options:
    -k=<k>      Length of sequences chosen for members of clique, computed each time as min(k, len(fasta_sequence)) [default: 31]
    -n=<n>      Number of sequences chosen for members of clique [default: 10]
    -B=<B>      Run bulk command with input file with paths to every clique to choose sequences for [default: 0]
    -h --help   Show this screen.
"""

from docopt import docopt
import pandas as pd
import numpy as np

# Reads a path to .tsv file, returns the bitmatrix, 
# indexes pointing to fasta files (sort of redundant)
# and clique member names in order relative to bitmatrix
def read_bitmatrix(path):
    tsv_file = pd.read_csv(path, sep='\t').sort_values('query_name')
    bitmat = tsv_file.values[:,1:]
    index = tsv_file.values[:,0]
    members = list(tsv_file.columns[1:])
    return(bitmat,index,members)

# From a bitmatrix, get a submatrix only consisting of unique
# vector forms founds in the original matrix
# Returns a reduced matrix consisting of every variation of 
# row vector, i.e. every way a node in CCDBG can be colored 
# differently, as well as a reference list which maps the 
# row index of a vector in the reduced matrix to the indices of
# that vector form's corresponding fasta sequences.
# Also returns members, a list of member names in order
# of appearance in bitmatrix columns, and finally
# a mapping from member indices to indices of sequences which
# they contain.
def reduce_bit_matrix(path):
    bitmat, index, members = read_bitmatrix(path)
    uniques = set() #set for our 'possibility-matrix'
    fasta_mapping = dict() #maps vector forms to every fasta sequence having that form
    max_available = 0
    for i in range(bitmat.shape[0]):
        if (bitmat[i] == 0).any(): #Weeds out vectors of pure 1's
            form = tuple(bitmat[i]) #Tuple for set functionality
            uniques.add(form) #Gather our vector forms into set
            if form in fasta_mapping:
                fasta_mapping[form].append(index[i]) #Add the fasta sequence's index to fasta_mapping for that vector form
            else:
                fasta_mapping[form] = [index[i]]
            max_available += 1
    N = len(uniques)
    M = len(members)
    reduced = np.zeros((N,M)) #Our new 'possibility-matrix'
    connected_to = dict() #Maps from index of member to indices of sequences in reduced possibility matrix appearing in that member
    ref_to_fasta = list() #We're gonna store our fasta sequence indexes in this list of lists, with index in the list matching to index in possibility matrix
    #!!! Be sure this sorting is okay, with respect to the order of relative indices !!!
    uniques = sorted(uniques, key=sum) #Sort our set based on degree of vector forms, making it into a list
    for i,form in enumerate(uniques): #For all unique vector forms in our bitmatrix (all vectors in new possibility matrix)
            reduced[i] = np.array(form) #Move it into an index of the possibility matrix
            ref_to_fasta.append(fasta_mapping[form]) #And log indices of all sequences taking that form in the same index in this
            for mem in np.where(reduced[i] == 1)[0]: #Find out which member(s) this vector form relates to
                if mem in connected_to:
                    connected_to[mem].append(i) #Log these vector forms as connected to those members
                else:
                    connected_to[mem] = [i]
    return(reduced, ref_to_fasta, members, connected_to, max_available)

def assign_sequences(reduced, fasta_indices, connected_to, iter_, available):
    M = reduced.shape[1]

    deg_1_span = np.where(np.sum(reduced, axis=1) == 2)[0][0]

    connections = np.zeros((M,M)) #Create connection matrix, logging how often each individual has been connected to another
    chosen = dict() #This will store our final selection of sequences and how they relate to members of clique
    for seq in range(deg_1_span): #Iterate over those having degree 1
        respective_member = int(np.where(reduced[seq] == 1)[0]) #Get the corresponding member related to vector
        for i, fasta in enumerate(fasta_indices[seq]): #Get the sub-list of every fasta sequence having this form
            if respective_member in chosen:
                chosen[respective_member].append(fasta) #Add all those fasta sequences to our chosen set
            else:
                chosen[respective_member] = [fasta]
            available -= 1
            connections[respective_member] += reduced[seq] #Log that we've found a sequence relating to that member
            i+=1
            if i == iter_: #If we've already reached our iter threshold (default 10), that's all we need...
                connected_to[respective_member].remove(seq)
                break

    while (sum(connections) < iter_).any() and available:
        focus_mem = np.argmin(np.sum(connections,axis=0)) #The member with fewest associated sequences so far...
        options = connected_to[focus_mem] #Get the index of the sequences associated with that member
        #Multiply the submatrix of reduced containing only associated sequences of member to the row vector relating to that
        #individual in the connections matrix. The lowest value in the resulting vector points to our best sequence to pick.
        #Because reduced is sorted, if many sequences are equally optimal, we will get lowest degree sequence of those
        chosen_index = np.argmin(np.matmul(reduced[options],connections[focus_mem])) 
        chosen_vec = options[chosen_index] #Get the index of that best sequence to pick
        chosen_fasta = fasta_indices[chosen_vec][0] #Pick any (the first) fasta sequence that has that vector form, get its index
        chosen.append(chosen_fasta) #Add that sequence index to our chosen sequence indices
        fasta_indices[chosen_vec].remove(chosen_fasta) #And remove that specific sequence from our pot of options
        if len(fasta_indices[chosen_vec]) == 0:
            for mem in np.where(reduced[chosen_vec] == 1)[0]:
                connected_to[mem].remove(chosen_vec)
        respective_members = np.where(reduced[chosen_vec] == 1)[0] #Find the members associated with that sequence...
        for respective_member in respective_members: #For each one...
            if respective_member in chosen:
                chosen[respective_member].append(chosen_fasta) #Log the sequence as associated with that member
            else:
                chosen[respective_member] = [chosen_fasta]
            available -= 1
            connections[respective_member] += reduced[chosen_vec] #And log new connections in connection matrix for every member associated
        
    return(chosen)
    

def parse_fasta(fasta_path):
        seqs = []
        s = []
        for line in open(fasta_path, 'r'):
            if line[0] == '>':
                if s:
                    seqs.append(''.join(s))
                s = []
            else:
                s.append(line.strip())
        seqs.append(''.join(s))
        return {'file': fasta_path, 'seqs': seqs}

#Writes results to a .seqs file, containing the mapping 
def print_to_file(chosen, mem, clique_path, seq_length):
    fastas = parse_fasta(clique_path+'.fasta')['seqs'] #Load our fasta sequences, with indices matching the indices used previously
    clique_id = clique_path.split('.')[0] #Create path to file with same prefix as other files
    with open(clique_id+'.seqs', 'w') as f:
        for mem_ind in chosen:
            f.write('>')
            f.write(mem[mem_ind].split('/')[-1].split('_genomic')[0]) #Get just NCBI name
            f.write('\n')
            for fas_ind in chosen[mem_ind]:
                fasta_seq = fastas[fas_ind][:min(seq_length, len(fastas[fas_ind]))] #Write sequences with at most given seq_length
                f.write(fasta_seq)
                f.write('\n')
            f.write('\n')

def choose_sequences(clique_path, k_mer_length, iter_):
    
    red, fas, mem, connected_to, max_available = reduce_bit_matrix(clique_path+'.tsv') #Get smaller 'possibility-matrix' and some connections
    chosen = assign_sequences(red, fas, connected_to, iter_, max_available) #Assigning indices of sequences to members
    
    print_to_file(chosen, mem, clique_path+'.fasta', k_mer_length) #Print results to file
    return

def load_paths(ref_path):
    with open(ref_path.strip()) as f:
        lines = f.readlines()
        return(lines)

def main(args):
    input_file = args['<input_file>']
    k_mer_length = int(args['-k'])
    iter_ = int(args['-n'])
    bulk = int(args['-B'])
    if bulk:
        cliques = load_paths(input_file)
        for clique in cliques:
            choose_sequences(clique, k_mer_length, iter_)
            print(clique, "done!")
    else:
        choose_sequences(input_file, k_mer_length, iter_)
        print(input_file, "done!")
    return(0)

if __name__ == "__main__":
    args = docopt(__doc__)
    main(args)