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

uint64_t convertToNumber(std::vector<int> colvec) {
    uint64_t sum = 0;
    int j = 0;
    for(int i = colvec.size() - 1; i >= 0; i--)
    {
        sum += colvec[i]<<j;
        j++;
    }
    return(sum);
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
Matrix readAndReduceBitMatrix(
    std::vector<std::string>* members,
    std::string clique_path, 
    khash_t(vector_u64)* fasta_mapping, 
    khash_t(vector_u64)* connected_to,
    std::vector<int>* goal_vector)
{
    int ret;
    khiter_t k;
    set<std::vector<int>> uniqueSet;

    std::ifstream tsv_data(clique_path+".tsv");
    std::string line;
    uint64_t fasta_ind;
    uint64_t vector_form;

    std::getline(tsv_data,line);
    std::vector<std::string> mems = splitString(line, '\t');
    mems = std::vector<std::string>(mems.begin() + 1, mems.end());
    // Ugly, temporary code
    for(auto member_name : mems) {
        std::size_t PosOfDir = member_name.find_last_of("/");
        std::size_t PosOfIDEnd = member_name.find_last_of("_");

        std::string NCBI_id = member_name.substr(PosOfDir+1, PosOfIDEnd);
        (*members).push_back(NCBI_id); // Add last substring to members
    }

    while(std::getline(tsv_data,line))
    {
        std::vector<int> line_split;
        std::vector<std::string> string_vector = splitString(line, '\t'); // Read lines in .tsv bitmatrix
        for(auto entry : string_vector) {
            line_split.push_back(stoi(entry));
        }
        fasta_ind = line_split[0]; // Get index of this fasta sequence
        line_split = std::vector<int>(line_split.begin() + 1, line_split.end()); // Split rest up into bit vector
        
        if(accumulate(line_split.begin(),line_split.end(),0) == line_split.size()) {
            continue; // We don't one pure 1's
        }
        addVectors(goal_vector, &line_split);
        uniqueSet.insert(line_split); // Insert our line into set of uniques
        vector_form = convertToNumber(line_split); // Transform bit vector to int as vector form identifier

        k = kh_get(vector_u64, fasta_mapping, vector_form); // Load fasta_mapping hash table

        if (k == kh_end(fasta_mapping))
        {
            k = kh_put(vector_u64, fasta_mapping, vector_form, &ret);
            kh_value(fasta_mapping, k) = new std::vector<uint64_t>;
        }
        kh_value(fasta_mapping, k)->push_back(fasta_ind); // Log this fasta sequence as belonging to vector form identifier
        
    }
    Matrix uniques;
    uniques.reserve(uniqueSet.size());
    for(auto vec : uniqueSet) {
        uniques.push_back(vec);
    }
    // Taken from https://stackoverflow.com/a/14419565
    std::sort(uniques.begin(), uniques.end(),
          [](const std::vector<int>& a, const std::vector<int>& b) {
        return accumulate(a.begin(), a.end(), 0) < accumulate(b.begin(), b.end(), 0);
    });
    for(int j = 0; j < uniques.size(); j++) { // For each of our unique vectors
        auto bit_vector = uniques[j];
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
uint64_t findFocusRow(int N, Matrix* uniques, std::vector<uint64_t>* option_rows, std::vector<int> connection_column) {
    size_t min_val = connection_column.size()*N;
    uint64_t focus_row = 0;
    for(auto opt : *option_rows) {
        auto this_val = dotProd((*uniques)[opt],connection_column);
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
    khash_t(vector_u64)* fasta_mapping,
    khash_t(vector_u64)* connected_to,
    Matrix* connections,
    Matrix* uniques)
{
    int ret;
    khiter_t k;
    auto chosen_vector = (*uniques)[chosen_vector_ind];
    int vector_form = convertToNumber(chosen_vector); // Get vector id number for this form
    k = kh_get(vector_u64, fasta_mapping, vector_form); // Find fasta sequences associated with that form
    auto value = kh_value(fasta_mapping, k); // Get the vector of those sequences
    uint64_t chosen_sequence = (*value)[0]; // Pick the first one (doesn't matter which)
    for(int i = 0; i < chosen_vector.size(); i++) {
        if(chosen_vector[i]){ // For every member taking that
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
                auto value = kh_value(connected_to, h); // Get a pointer to its value
                std::remove((*value).begin(), (*value).end(), chosen_vector_ind); // Remove connection from member to vector from (ONLY for singles)
            }
        }
    }
    std::remove((*value).begin(), (*value).end(), chosen_sequence);
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
    khash_t(vector_u64)* fasta_mapping, 
    khash_t(vector_u64)* connected_to,
    khash_t(vector_u64)* chosen,
    Matrix* connections,
    Matrix* uniques,
    std::vector<int>* goal_vector)
{
    int ret;
    khiter_t k;
    Matrix::iterator itr;
    int vector_form;
    int keep_track = 0;
    // First pass with degree 1
    for (itr = (*uniques).begin(); itr != (*uniques).end(); itr++)
    {
        if(accumulate((*itr).begin(), (*itr).end(), 0) > 1) {
            break; // Now we've reached higher degree vectors in uniques (matrix sorted by degree)
        }
        for(size_t i = 0; i < (*itr).size(); i++) {
            if((*itr)[i] && (*connections)[i][i] < N+1) { // For 1 in results
                vector_form = convertToNumber((*itr)); // Get its corresponding vector form
                k = kh_get(vector_u64, fasta_mapping, vector_form); // Find fasta sequences associated with that form
                auto value = kh_value(fasta_mapping, k); // Get the vector of those sequences
                k = kh_get(vector_u64, chosen, i); // Find the member's location in chosen hash table
                if (k == kh_end(chosen))
                {
                    k = kh_put(vector_u64, chosen, i, &ret);
                    kh_value(chosen, k) = new std::vector<uint64_t>;
                }
                int count = 0;
                for(auto val : *value) { // For indices of fasta sequences with that vector form
                    kh_value(chosen, k)->push_back(val); // Add them to the chosen hash table for that member
                    count++;
                    (*connections)[i][i]++; // Log it to connection matrix
                    if(count == N) { // We want max N sequences
                        break;
                    }
                }
                k = kh_get(vector_u64, connected_to, i); // Then find the location of the member in connected_to hash table
                auto value_connected = kh_value(connected_to, k); // Get a pointer to its value
                std::remove((*value_connected).begin(), (*value_connected).end(), keep_track); // Remove connection from member to vector from (ONLY for singles)
                break;
            }
        }
        keep_track++;
    }
    while(goalsRemaining(connections, goal_vector, N)) {
        // Second pass for higher degrees
        int focus_member = findFocusMember(connections, goal_vector, N); // Find member with lowest number of connections
        k = kh_get(vector_u64, connected_to, focus_member); // Then find the location of the member in connected_to hash table
        auto options = kh_value(connected_to, k); // Get a pointer to its value, showing its associated vector forms
        uint64_t focus_row = findFocusRow(N, uniques, options, (*connections)[focus_member]);
        logSequenceToMembers(focus_row, chosen, fasta_mapping, connected_to, connections, uniques);
    }
    return(chosen);
}

std::vector<pair<int, Kmer>> parseFastaSeq(size_t k, std::string fasta_file, size_t size) {

    Kmer::set_k(k);
    std::vector<pair<int, Kmer>> fasta_sequences;
    fasta_sequences.reserve(size);
    // Taken from https://stackoverflow.com/a/35251990
    // Will be substituted by bifrost FASTX_Parser once I know it works.
    std::string line, id, DNA_sequence;
    Kmer seq;
    std::ifstream infile(fasta_file);
    std::string filename;

    if (!infile.good()) {
        std::cerr << "Error opening: " << fasta_file << ". Path is incorrect." << std::endl;
        return(fasta_sequences);
    }
    // Don't loop on good(), it doesn't allow for EOF!!

    //    while (std::getline(input, line).good()) {
    while (std::getline(infile, line)) {

        // line may be empty so you *must* ignore blank lines
        // or you have a crash waiting to happen with line[0]
        if(line.empty())
            continue;

        if (line[0] == '>') {
            // output previous line before overwriting id
            // but ONLY if id actually contains something
            if(!id.empty()) {
                DNA_sequence = DNA_sequence.substr(0,min(k, DNA_sequence.size()));
                seq = Kmer(DNA_sequence.c_str()).rep(); // rep() gets canonical Kmer
                fasta_sequences.push_back(std::make_pair(stoi(id), seq));
            }
            id = line.substr(1);
            DNA_sequence.clear();
        }
        else {
            DNA_sequence += line;
        }
    }
    // This is the default, but it's put here to be clear.
    std::sort(fasta_sequences.begin(), fasta_sequences.end(),
          [](const pair<int, Kmer>& a, const pair<int, Kmer>& b) {
        return a.first < b.first;
    });
    return(fasta_sequences);
}

void writeFileToSeq(
    int k,
    std::string infile,
    size_t clique_size_unitigs, 
    khash_t(vector_u64)* chosen,
    std::vector<std::string>* members) 
{
    auto fasta_sequences = parseFastaSeq(k, infile+".fasta", clique_size_unitigs);
    
    std::size_t PosOfFileExt = infile.find_last_of(".");
    std::string clique_id = infile.substr(0, PosOfFileExt);
    std::ofstream seqfile(clique_id + ".minced");

    for (khiter_t k = kh_begin(chosen); k != kh_end(chosen); ++k) {
        if (kh_exist(chosen, k)) {
            auto key = kh_key(chosen, k);
            auto val = kh_val(chosen, k);

            if (val->size() > 1)
            {
                int i = 0;
                for(auto seq_index : *val) {
                    seqfile << '>' << (*members)[key] << ';' << infile << '!' << '\n';
                    seqfile << fasta_sequences[seq_index].second.toString() << '\n';
                    i++;
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
    std::ofstream seqfile(clique_id + ".minced");
    seqfile << "None\n" << '\n' << "Clique " << infile << " has no distinguishing kmers to be found";
    seqfile.close();
}

void runSequenceFind(int k, int N, size_t clique_size_members, size_t clique_size_unitigs, std::string infile) {
    Matrix connections;
    for(int i = 0; i < clique_size_members; i++) {
        std::vector<int> vector1(clique_size_members, 0);
        connections.push_back(vector1);
    }

    khash_t(vector_u64)* fasta_mapping = kh_init(vector_u64);
    khash_t(vector_u64)* connected_to = kh_init(vector_u64);
    khash_t(vector_u64)* chosen = kh_init(vector_u64);
    std::vector<int> goal_vector(clique_size_members, 0);
    std::vector<std::string> members;

    Matrix uniques = readAndReduceBitMatrix(&members, infile, fasta_mapping, connected_to, &goal_vector);
    if(uniques.size() == 0) {
        writeFailedCliqueToFile(infile, &members);
        std::cout << "Clique " << infile << " has no distinguishing kmers to be found" << std::endl;
        return;
    }
    
    chosen = chooseSequences(N, fasta_mapping, connected_to, chosen, &connections, &uniques, &goal_vector);
    writeFileToSeq(k, infile, clique_size_unitigs, chosen, &members);
    delete fasta_mapping;
    delete connected_to;
    delete chosen;
}