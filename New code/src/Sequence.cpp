#include <fstream>
#include <getopt.h>
#include <vector>
#include <set>
#include <numeric>
#include <map>
#include <iostream>
#include <bifrost/CompactedDBG.hpp>
#include "khash.h"

KHASH_MAP_INIT_INT64(vector_u64, std::vector<uint64_t>*);

auto cmp = [](std::vector<int> a, std::vector<int> b) 
{
    return accumulate(a.begin(),a.end(),0) - accumulate(b.begin(),b.end(),0);
};

using UniqueSet = std::set<std::vector<int>, decltype(cmp)>;
using Matrix = std::vector<std::vector<int>>;

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.txt> | <in.clique>) ...\n\n"
        "Options:\n"
        
        "  -sw         Size of clique, i.e. number of genomes.\n"
        "  -d <path>   Destination directory for sequence files.\n"
        "  -k          Size of sequences [default: 31]\n"
        "  -N          Number of sequences per genome in clique [default: 10]";
    printf(s, name);
}

auto splitString(std::string_view in, char sep) {
    std::vector<int> r;
    r.reserve(std::count(in.begin(), in.end(), sep) + 1); // optional
    for (auto p = in.begin();; ++p) {
        auto q = p;
        p = std::find(p, in.end(), sep);
        r.emplace_back(q, p);
        if (p == in.end())
            return r;
    }
}

auto splitMemberString(std::string_view in, char sep) {
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

// Height from opt.filename_ref_in.size() and width from dbg.size() in makeBitMatrix.cpp 
Matrix* readAndReduceBitMatrix(
    std::vector<std::string>* members,
    std::string clique_path, 
    khash_t(vector_u64)* fasta_mapping, 
    khash_t(vector_u64)* connected_to,
    std::vector<int>* goal_vector)
{
    int ret;
    khiter_t k;
    UniqueSet uniqueSet(cmp);

    std::fstream tsv_data(clique_path+".tsv", std::ios::in);
    std::istream_iterator<uint64_t> start(tsv_data), end;

    std::string line;
    std::vector<int> line_split;
    uint64_t fasta_ind;
    uint64_t vector_form;

    tsv_data >> line; // Getting rid of the first line
    std::vector<std::string> mems = splitMemberString(line, '\t');
    mems = std::vector<std::string>(mems.begin() + 1, mems.end());
    // Ugly, temporary code
    for(auto member_name : mems) {
        istringstream f(member_name);
        string s;    
        while (getline(f, s, '/')) { // Iterating over substrings delimitered by '/'
        }
        (*members).push_back(s); // Add last substring to members
    }
    uint64_t i = 0;

    while(std::getline(tsv_data,line))
    {
        line_split = splitString(line, '\t'); // Read lines in .tsv bitmatrix
        fasta_ind = line_split[0]; // Get index of this fasta sequence
        line_split = std::vector<int>(line_split.begin() + 1, line_split.end()); // Split rest up into bit vector
        if(accumulate(line_split.begin(),line_split.end(),0) == line_split.size()) {
            continue; // We don't one pure 1's
        }
        (*goal_vector) = addVectors((*goal_vector), line_split);
        uniqueSet.insert(line_split); // Insert our line into set of uniques
        vector_form = convertToNumber(line_split); // Transform bit vector to int as vector form identifier

        k = kh_get(vector_u64, fasta_mapping, vector_form); // Load fasta_mapping hash table

        if (k == kh_end(fasta_mapping))
        {
            k = kh_put(vector_u64, fasta_mapping, vector_form, &ret);
            kh_value(fasta_mapping, k) = new std::vector<uint64_t>;
        }
        kh_value(fasta_mapping, k)->push_back(fasta_ind); // Log this fasta sequence as belonging to vector form identifier
        

        for(int i = 0; i < line_split.size(); i++) // Iterate through our bit vector
        {
            if(line_split[i]) { // If i = 1 in vector...
                k = kh_get(vector_u64, connected_to, i);
                if (k == kh_end(connected_to))
                {
                    k = kh_put(vector_u64, connected_to, i, &ret);
                    kh_value(connected_to, k) = new std::vector<uint64_t>;
                }
                kh_value(connected_to, k)->push_back(vector_form); // Log this vector as being connected to member i
            }
        }
        
    }
    Matrix uniques;
    for(auto vec : uniqueSet) {
        uniques.push_back(vec);
    }
    return(&uniques);
}

// WAY TO COMBINE FINDFOCUSMEMBER AND GOALSREMAINING?

// Computes column_sum of the connections matrix to determine which member has fewest connections
size_t findFocusMember(Matrix* connections, std::vector<int>* goal_vector, int N) {
    int max_missing = 0;
    int focus_mem;
    for(int i = 0; i < (*connections).size(); i++) {
        int sum = 0;
        int goal = std::min((*goal_vector)[i], N);
        for(int j = 0; j < (*connections).size(); j++) {
            sum+=(*connections)[j][i];
        }
        if((goal - sum) > max_missing) {
            max_missing = sum;
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
uint64_t findFocusRow(Matrix* uniques, std::vector<uint64_t>* option_rows, std::vector<int> connection_column) {
    size_t min_val = connection_column.size();
    uint64_t focus_column;
    for(auto opt : *option_rows) {
        auto this_val = dotProd((*uniques)[opt],connection_column);
        if(this_val < min_val) {
            min_val = this_val;
            focus_column = opt;
        }
    }
    return(focus_column);
}

std::vector<int> addVectors(std::vector<int> vec1, std::vector<int> vec2) {
    if(vec1.size() != vec2.size()) {
        std::cout << "Error: vector sizes not compatible for: \n";
        std::cout << vec1.size() << vec2.size() << std::endl;
    }
    std::vector<int> res = vec1;
    for(int i = 0; i < vec2.size(); i++) {
        res[i] += vec2[i];
    }
    return(res);
}

void logSequenceToMembers(
    std::vector<int> chosen_vector,
    khash_t(vector_u64)* chosen,
    khash_t(vector_u64)* fasta_mapping,
    khash_t(vector_u64)* connected_to,
    Matrix* connections)
{
    int ret;
    khiter_t k;
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
            (*connections)[i] = addVectors((*connections)[i], chosen_vector);

            if((*value).size() == 1) { // If we have no more sequences associated with this vector form, except for this one
                khiter_t h;
                h = kh_get(vector_u64, connected_to, i); // Then find the location of the member in connected_to hash table
                auto value = kh_value(connected_to, h); // Get a pointer to its value
                std::remove((*value).begin(), (*value).end(), vector_form); // Remove connection from member to vector from (ONLY for singles)
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
    // First pass with degree 1
    for (itr = (*uniques).begin(); itr != (*uniques).end(); itr++)
    {
        if(accumulate((*itr).begin(), (*itr).end(), 0) > 1) {
            break;
        }
        for(size_t i = 0; i < (*itr).size(); i++) {
            if((*itr)[i] && (*connections)[i][i] < N+1) { // For 1 in results
                (*connections)[i][i]++; // Log it to connection matrix
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
                    if(count == N) { // We want max N sequences
                        break;
                    }
                }
                k = kh_get(vector_u64, connected_to, i); // Then find the location of the member in connected_to hash table
                auto value = kh_value(connected_to, k); // Get a pointer to its value
                std::remove((*value).begin(), (*value).end(), vector_form); // Remove connection from member to vector from (ONLY for singles)
                break;
            }
        }
    }
    while(goalsRemaining(connections, goal_vector, N)) {
        // Second pass for higher degrees
        int focus_member = findFocusMember(connections, goal_vector, N); // Find member with lowest number of connections
        k = kh_get(vector_u64, connected_to, focus_member); // Then find the location of the member in connected_to hash table
        auto options = kh_value(connected_to, k); // Get a pointer to its value, showing its associated vector forms
        uint64_t focus_row = findFocusRow(uniques, options, (*connections)[focus_member]);
        logSequenceToMembers((*uniques)[focus_row], chosen, fasta_mapping, connected_to, connections);
    }
    return(chosen);
}

std::vector<std::string> parseFastaSeq(size_t k, std::string fasta_file, size_t size) {
    std::vector<std::string> fasta_sequences;
    fasta_sequences.reserve(size);
    // Taken from https://stackoverflow.com/a/35251990
    // Will be substituted by bifrost FASTX_Parser once I know it works.
    std::string line, id, DNA_sequence;
    std::ifstream infile(fasta_file);
    std::string filename;

    if (!infile.good()) {
        std::cerr << "Error opening: " << fasta_file << " . You have failed." << std::endl;
        return;
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
            if(!id.empty())
                fasta_sequences.push_back(DNA_sequence.substr(0,min(k, DNA_sequence.size())));

            id = line.substr(1);
            DNA_sequence.clear();
        }
        else {//  if (line[0] != '>'){ // not needed because implicit
            DNA_sequence += line;
        }
    }
    return(fasta_sequences);
}

int main(int argc, char* argv[]) 
{
    if (argc == 1)
    {
        print_usage(argv[0]);
        exit(1);
    }

    std::string infile = argv[argc-1];
    std::string dirpath = "";
    uint16_t k = 31;
    uint16_t N = 10;
    uint16_t clique_size;
    size_t size_of_sequence_file;

    int opt;
    while ((opt = getopt(argc, argv, "s:d:k:N:")) != -1)
    {
        switch (opt)
        {
            case 's':
                clique_size = atoi(optarg);
                break;
            case 'd':
                dirpath = optarg;
                if (dirpath.back() != '/')
                    dirpath += '/';
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'N':
                N = atoi(optarg);
                break;
        }
    }
    
    Matrix connections;
    for(int i = 0; i < clique_size; i++) {
        std::vector<int> vector1(clique_size, 0);
        connections.push_back(vector1);
    }

    khash_t(vector_u64)* fasta_mapping = kh_init(vector_u64);
    khash_t(vector_u64)* connected_to = kh_init(vector_u64);
    khash_t(vector_u64)* chosen = kh_init(vector_u64);
    std::vector<int> goal_vector(clique_size, 0);
    std::vector<std::string> members;

    Matrix* uniques = readAndReduceBitMatrix(&members, infile, fasta_mapping, connected_to, &goal_vector);
    chosen = chooseSequences(N, fasta_mapping, connected_to, chosen, &connections, uniques, &goal_vector);

    // Switch this out for bifrost FASTX_Parser when I've got things at least running...
    auto fasta_sequences = parseFastaSeq(k, infile+".fasta", size_of_sequence_file);
    std::ofstream seqfile(infile + ".seq");

    for (khiter_t k = kh_begin(chosen); k != kh_end(chosen); ++k) {
        if (kh_exist(chosen, k)) {
            auto key = kh_key(chosen, k);
            auto val = kh_val(chosen, k);

            if (val->size() > 1)
            {
                seqfile << '>' << members[key] << '\n';
                for(auto seq_index : *val) {
                    seqfile << fasta_sequences[seq_index] << '\n';
                }
                seqfile << '\n';
            }
        }
    }
    seqfile.close();
}