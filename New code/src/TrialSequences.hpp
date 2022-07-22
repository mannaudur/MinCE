#include <fstream>
#include <vector>
#include <set>
#include <numeric>
#include <map>
#include <iostream>
#include <bifrost/ColoredCDBG.hpp>
#include <bifrost/CompactedDBG.hpp>
#include "khash.h"

KHASH_MAP_INIT_INT64(vector_u64, std::vector<uint64_t>*);

using Matrix = std::vector<std::vector<int>>;
using Uniques = std::vector<std::pair<std::vector<int>,std::vector<uint64_t>>>;

std::vector<std::string> splitString(std::string in, char sep) {
    std::vector<std::string> r;
    r.reserve(std::count(in.begin(), in.end(), sep) + 1); // optional
    for (auto p = in.begin();; ++p) {
        auto q = p;
        p = std::find(p, in.end(), sep);
        r.emplace_back(q, p);
        if (p == in.end())
            return r;
    }
}

std::vector<int> splitInt(std::string in, char sep) {
    std::vector<std::string> r;
    std::vector<int> line_split;
    r.reserve(std::count(in.begin(), in.end(), sep) + 1); // optional
    for (auto p = in.begin();; ++p) {
        auto q = p;
        p = std::find(p, in.end(), sep);
        r.emplace_back(q, p);
        if (p == in.end()) {
            for(auto entry : r) {
                line_split.push_back(stoi(entry));
            }
            return line_split;
        }
    }
}

void addVectors(std::vector<int>* vec1, std::vector<int>* vec2) {
    if((*vec1).size() != (*vec2).size()) {
        std::cout << "Error: vector sizes not compatible for: \n";
        std::cout << (*vec1).size() << (*vec2).size() << std::endl;
    }
    for(int i = 0; i < (*vec2).size(); i++) {
        (*vec1)[i] += (*vec2)[i];
    }
}

// Height from opt.filename_ref_in.size() and width from dbg.size() in Distinguish.cpp 
Uniques readAndReduceBitMatrix(
    std::vector<std::string>* members,
    std::string clique_path,
    khash_t(vector_u64)* connected_to,
    std::vector<int>* goal_vector)
{
    int ret;
    khiter_t k;

    std::ifstream tsv_data(clique_path+".tsv");
    std::string line;
    uint64_t fasta_ind;
    std::vector<int> bit_vector;

    std::getline(tsv_data,line);
    std::vector<std::string> mems = splitString(line, '\t');
    mems = std::vector<std::string>(mems.begin() + 1, mems.end());
    
    for(auto member_name : mems) {
        std::size_t IDEnd = member_name.find("_genomic");
        std::size_t IDBegin = member_name.find_last_of('/');
        std::string genome_NCBI = member_name.substr(0, IDEnd);
        genome_NCBI = genome_NCBI.substr(IDBegin+1);

        (*members).push_back(genome_NCBI); // Add substring having only NCBI identifier to members
    }

    Uniques uniques;

    while(std::getline(tsv_data,line))
    {
        std::vector<int> line_split = splitInt(line, '\t'); // Read lines in .tsv bitmatrix
        fasta_ind = line_split[0]; // Get index of this fasta sequence
        bit_vector = std::vector<int>(line_split.begin() + 1, line_split.end()); // Split rest up into bit vector
        
        if(accumulate(bit_vector.begin(),bit_vector.end(),0) == bit_vector.size()) {
            continue; // We don't one pure 1's
        }
        bool found = false;
        addVectors(goal_vector, &bit_vector);
        for(int i = 0; i < uniques.size(); i++) {
            if(uniques[i].first == bit_vector) {
                uniques[i].second.push_back(fasta_ind);
                found = true;
            }
        }
        if(!found) {
            std::vector<uint64_t> fasta_vector;
            fasta_vector.push_back(fasta_ind);
            uniques.push_back(std::make_pair(bit_vector, fasta_vector)); // Insert our line into set of uniques
        }
        
    }

    // Taken from https://stackoverflow.com/a/14419565
    std::sort(uniques.begin(), uniques.end(),
          [](const std::pair<std::vector<int>, std::vector<uint64_t>>& a, 
          const std::pair<std::vector<int>, std::vector<uint64_t>>& b) {
        return accumulate(a.first.begin(), a.first.end(), 0) < accumulate(b.first.begin(), b.first.end(), 0);
    });
    for(int j = 0; j < uniques.size(); j++) { // For each of our unique vectors
        auto bit_vector = uniques[j].first;
        for(int i = 0; i < bit_vector.size(); i++) // Iterate through our bit vector
        {
            if(bit_vector[i]) { // If i = 1 in vector...
                k = kh_get(vector_u64, connected_to, i);
                if (k == kh_end(connected_to))
                {
                    k = kh_put(vector_u64, connected_to, i, &ret);
                    kh_value(connected_to, k) = new std::vector<uint64_t>;
                }
                kh_value(connected_to, k)->push_back(j); // Log this vector as being connected to member i
            }
        }
    }

    return(uniques);
}

// Computes column_sum of the connections matrix to determine which member has fewest connections
size_t findFocusMember(Matrix* connections, std::vector<int>* goal_vector, int N) {
    int max_missing = 0;
    int focus_mem = 0;
    for(int i = 0; i < (*connections).size(); i++) {
        int sum = 0;
        int goal = std::min((*goal_vector)[i], N);
        for(auto num : (*connections)[i]) {
            sum+=num;
        }
        if((goal - sum) > max_missing) {
            max_missing = (goal - sum);
            focus_mem = i;
        }
    }
    return(focus_mem);
}

// Simple dot product calculation
int dotProd(std::vector<int> unique_row, std::vector<int> connection_column) {
    int sum = 0;
    for(int i = 0; i < unique_row.size(); i++) {
        sum += (unique_row[i]*connection_column[i]);
    }
    return(sum);
}

// Computes dot product of chosen focus column of connections and every possible vector row in possibility matrix
// As possibility matrix, uniques, is sorted by rowsum, we will get the lowest degree vector in the case of a tie
// for lowest rowsum.
uint64_t findFocusRow(int N, Uniques* uniques, std::vector<uint64_t>* option_rows, std::vector<int> connection_column) {
    size_t min_val = connection_column.size()*N;
    uint64_t focus_row = 0;
    for(auto opt : *option_rows) {
        auto this_val = dotProd((*uniques)[opt].first,connection_column);
        if(this_val < min_val) {
            min_val = this_val;
            focus_row = opt;
        }
    }
    return(focus_row);
}


void logSequenceToMembers(
    uint64_t chosen_vector_ind,
    khash_t(vector_u64)* chosen,
    khash_t(vector_u64)* connected_to,
    Matrix* connections,
    Uniques* uniques,
    size_t iter_)
{
    int ret;
    khiter_t k;
    auto chosen_vector = (*uniques)[chosen_vector_ind].first;
    auto value = &(*uniques)[chosen_vector_ind].second; // Get the vector of those sequence ids
    uint64_t chosen_sequence; 
    size_t size_of_options = value->size();
    for(int j = 0; j < min(iter_, size_of_options); j++) {
        chosen_sequence = value->back(); // Choose the last sequence in the vector
        for(int i = 0; i < chosen_vector.size(); i++) {
            if(chosen_vector[i]){ // For every member taking that sequence
                k = kh_get(vector_u64, chosen, i); // Find the member's location in chosen hash table
                if (k == kh_end(chosen))
                {
                    k = kh_put(vector_u64, chosen, i, &ret);
                    kh_value(chosen, k) = new std::vector<uint64_t>;
                }
                kh_value(chosen, k)->push_back(chosen_sequence);
                
                addVectors(&((*connections)[i]), &chosen_vector);

                if(value -> size() == 1) { // If we have no more sequences associated with this vector form, except for this one
                    khiter_t h;
                    h = kh_get(vector_u64, connected_to, i); // Then find the location of the member in connected_to hash table
                    auto value_connected = kh_value(connected_to, h); // Get a pointer to its value
                    std::vector<uint64_t>::iterator it;
                    it = std::find(value_connected->begin(), value_connected->end(), chosen_vector_ind);
                    value_connected->erase((*value_connected).begin() + (it-value_connected->begin())); // Remove connection from member to vector from (ONLY for singles)
                }
            }
        }
        (*uniques)[chosen_vector_ind].second.pop_back();
    }
}

bool goalsRemaining(Matrix* connections, std::vector<int>* goal_vector, int N) {
    for(int i = 0; i < (*connections).size(); i++) {
            int sum = 0;
            int goal = std::min((*goal_vector)[i], N);
            for(int j = 0; j < (*connections).size(); j++) {
                sum+=(*connections)[j][i];
            }
            if((goal - sum) > 0) {
                return(true);
            }
        }
    return(false);
}

khash_t(vector_u64)* chooseSequences(
    int N,
    khash_t(vector_u64)* connected_to,
    khash_t(vector_u64)* chosen,
    Matrix* connections,
    Uniques* uniques,
    std::vector<int>* goal_vector)
{
    int ret;
    khiter_t k;
    // First pass with degree 1
    for (int row_num = 0; row_num < uniques->size(); row_num++)
    {
        if(accumulate((*uniques)[row_num].first.begin(), (*uniques)[row_num].first.end(), 0) > 1) {
            break; // Now we've reached higher degree vectors in uniques (matrix sorted by degree)
        }
        auto row_vector = (*uniques)[row_num];

        for(size_t i = 0; i < row_vector.first.size(); i++) {
            if(row_vector.first[i] && (*connections)[i][i] < N+1) { // For 1 in results
                auto value = &row_vector.second; // Get the vector of those sequences
                k = kh_get(vector_u64, chosen, i); // Find the member's location in chosen hash table
                if (k == kh_end(chosen))
                {
                    k = kh_put(vector_u64, chosen, i, &ret);
                    kh_value(chosen, k) = new std::vector<uint64_t>;
                }

                while(value->size()) { // For indices of fasta sequences with that vector form
                    kh_value(chosen, k)->push_back(value->back()); // Add them to the chosen hash table for that member
                    value->pop_back();
                    if(++(*connections)[i][i] > N-1) // Log it to connection matrix
                        break;
                }
                k = kh_get(vector_u64, connected_to, i); // Then find the location of the member in connected_to hash table
                auto value_connected = kh_value(connected_to, k); // Get a pointer to its value
                std::vector<uint64_t>::iterator it;
                it = std::find(value_connected->begin(), value_connected->end(), row_num);
                value_connected->erase((*value_connected).begin() + (it-value_connected->begin())); 
                break;
            }
        }
    }
    while(goalsRemaining(connections, goal_vector, N)) {
        // Second pass for higher degrees
        int focus_member = findFocusMember(connections, goal_vector, N); // Find member with lowest number of connections
        k = kh_get(vector_u64, connected_to, focus_member); // Then find the location of the member in connected_to hash table
        auto options = kh_value(connected_to, k); // Get a pointer to its value, showing its associated vector forms
        uint64_t focus_row = findFocusRow(N, uniques, options, (*connections)[focus_member]);
        logSequenceToMembers(focus_row, chosen, connected_to, connections, uniques, 1);
    }
    // Last pass to check whether any members have not been separated through the process above
    // For clarity: connections[i][i] is the diagonal, which counts how many times member i has been issued a sequence.
    // If any cell in the corresponding vector has the same value, those two members are not separated,
    // so we try to find at least 3 sequences to separate them, if we can.
    std::vector<uint64_t> final_adds;
    const size_t C = 3;
    for(int i = 0; i < connections->size()-1; i++) {
        for(int j = i+1; j < connections->size(); j++) {
            if((*connections)[i][i] == (*connections)[i][j]) 
            {
                k = kh_get(vector_u64, connected_to, i); 
                if (k != kh_end(connected_to))
                {
                    auto options = kh_value(connected_to, k);
                    for(auto opt : (*options)) {
                        if(!(*uniques)[opt].first[j])
                            final_adds.push_back(opt);
                    }
                    for(int M = 0; M < min(final_adds.size(), C); M++) {
                        uint64_t focus_row = findFocusRow(N, uniques, &final_adds, (*connections)[i]);
                        logSequenceToMembers(focus_row, chosen, connected_to, connections, uniques, 3);
                    }
                }
                final_adds.clear();

                k = kh_get(vector_u64, connected_to, j); 
                if (k != kh_end(connected_to))
                {
                    auto options = kh_value(connected_to, k);
                    for(auto opt : (*options)) {
                        if(!(*uniques)[opt].first[i])
                            final_adds.push_back(opt);
                    }
                    for(size_t M = 0; M < min(final_adds.size(), C); M++) {
                        uint64_t focus_row = findFocusRow(N, uniques, &final_adds, (*connections)[j]);
                        logSequenceToMembers(focus_row, chosen, connected_to, connections, uniques, 3);
                    }
                }
                final_adds.clear();
            }
        }
    }
    return(chosen);
}

std::vector<pair<int, Kmer>> parseFastaSeq(size_t k, std::string fasta_file, size_t size) 
{
    Kmer::set_k(k);
    std::vector<pair<int, Kmer>> fasta_sequences;
    fasta_sequences.reserve(size);
    std::string line, id, DNA_sequence;
    Kmer sequence;
    std::vector<std::string> vec;
    string seq;
    size_t file_id;
    vec.push_back(fasta_file);
    FileParser f(vec);
    uint i = 0;

    while(f.read(seq,file_id)) {
        DNA_sequence = seq.substr(0,min(k, seq.size()));
        sequence = Kmer(DNA_sequence.c_str()).rep(); // rep() gets canonical Kmer
        fasta_sequences.push_back(std::make_pair(i++, sequence));
    }
    // This is the default, but it's put here to be clear.
    std::sort(fasta_sequences.begin(), fasta_sequences.end(),
          [](const pair<int, Kmer>& a, const pair<int, Kmer>& b) {
        return a.first < b.first;
    });
    return(fasta_sequences);
}

void writeFileToSeq(
    int k_,
    std::string infile,
    size_t clique_size_unitigs, 
    khash_t(vector_u64)* chosen,
    std::vector<std::string>* members) 
{
    auto fasta_sequences = parseFastaSeq(k_, infile+".fasta", clique_size_unitigs);
    
    std::size_t PosOfFileExt = infile.find_last_of(".");
    std::string clique_id = infile.substr(0, PosOfFileExt);
    std::ofstream seqfile(clique_id + ".sequences");

    for (khiter_t k = kh_begin(chosen); k != kh_end(chosen); ++k) {
        if (kh_exist(chosen, k)) {
            auto key = kh_key(chosen, k);
            auto val = kh_val(chosen, k);

            if (val->size() > 1)
            {
                seqfile << '>' << (*members)[key] << '\n';
                for(auto seq_index : *val) {
                    seqfile << fasta_sequences[seq_index].second.toString() << '\n';
                }
                seqfile << '\n';
            }
        }
    }
    seqfile.close();
}

void writeFailedCliqueToFile(std::string infile, std::vector<std::string>* members) {
    std::size_t PosOfFileExt = infile.find_last_of(".");
    std::string clique_id = infile.substr(0, PosOfFileExt);
    std::ofstream seqfile(clique_id + ".sequences");
    for(auto mem : (*members)) {
        seqfile << '>' << mem << std::endl;
    }
    seqfile.close();
}

void runSequenceFind(int k, int N, size_t clique_size_members, size_t clique_size_unitigs, std::string infile) {
    Matrix connections;
    for(int i = 0; i < clique_size_members; i++) {
        std::vector<int> vector1(clique_size_members, 0);
        connections.push_back(vector1);
    }

    khash_t(vector_u64)* connected_to = kh_init(vector_u64);
    khash_t(vector_u64)* chosen = kh_init(vector_u64);
    std::vector<int> goal_vector(clique_size_members, 0);
    std::vector<std::string> members;
    Uniques uniques = readAndReduceBitMatrix(&members, infile, connected_to, &goal_vector);
    if(uniques.size() == 0) {
        writeFailedCliqueToFile(infile, &members);
        std::cout << "Clique " << infile << " has no distinguishing kmers to be found" << std::endl;
        return;
    }
    
    chosen = chooseSequences(N, connected_to, chosen, &connections, &uniques, &goal_vector);
    writeFileToSeq(k, infile, clique_size_unitigs, chosen, &members);
}