#include "khash.h"
#include <string>
#include <cstdint>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <bifrost/CompactedDBG.hpp>
#include <set>
#include <cstdint>

KHASH_MAP_INIT_INT64(vector_u64, std::vector<uint64_t>*);
using Index = std::vector<std::pair<std::string, int> >;

std::vector<std::string> readPathFile(std::string ifpath)
{
    std::vector<std::string> fnames;
    std::fstream fin;
    fin.open(ifpath, std::ios::in);

    if (fin.is_open())
    {
        std::string fname;
        while (getline(fin, fname))
        {
            fnames.push_back(fname);
        }
    }
    return fnames;
}

Index read_index(const char *index_filename, std::map<std::string, uint64_t> *index_mapping)
{
    int index_id;
    std::string genome;
    std::string clique;
    std::fstream fs(index_filename, std::ios::in);
    std::string line;
    Index index;
    while (std::getline(fs, line))
    {
        std::stringstream ss(line);
        ss >> index_id >> genome >> clique;

        std::size_t IDEnd = genome.find("_genomic");
        std::size_t IDBegin = genome.find_last_of('/');
        std::string genome_NCBI = genome.substr(0, IDEnd);
        genome_NCBI = genome_NCBI.substr(IDBegin+1);

        if (clique == "NULL")
            index.push_back(std::make_pair(genome_NCBI, -1));
        else {
            index.push_back(std::make_pair(genome_NCBI, 0));
            (*index_mapping)[genome_NCBI] = index_id;
        }
    }
    return index;
}

void addSeqToHashMap(std::string path, khash_t(vector_u64)* seq_locator, std::map<std::string, uint64_t> *index_mapping) {
    std::string line;
    std::string id;
    uint64_t MinCE_id;
    std::vector<uint64_t> hash_sequences;
    std::ifstream infile(path);
    khiter_t k;
    int ret;
    while (std::getline(infile, line)) {

        // line may be empty so you *must* ignore blank lines
        // or you have a crash waiting to happen with line[0]
        if(line.empty())
            continue;

        if (line[0] == '>') {
            // output previous line before overwriting id
            // but ONLY if id actually contains something
            if(!id.empty()) {
                for(auto hash : hash_sequences) {
                    k = kh_get(vector_u64, seq_locator, hash);

                    if (k == kh_end(seq_locator))
                    {
                        k = kh_put(vector_u64, seq_locator, hash, &ret);
                        kh_value(seq_locator, k) = new std::vector<uint64_t>;
                    }
                    kh_value(seq_locator, k)->push_back(MinCE_id);
                }
            }
            id = line.substr(1);
            MinCE_id = (*index_mapping)[id];
            hash_sequences.clear();
        }
        else 
        {
            Kmer seq = Kmer(line.c_str());
            hash_sequences.push_back(seq.hash());
        }
    }
}

void removeDuplicatesAndLogNoOfSeq(khash_t(vector_u64)* seq_locator, Index* index) {
    for (khiter_t k = kh_begin(seq_locator);
         k != kh_end(seq_locator);
         ++k)
    {
        if (kh_exist(seq_locator, k))
        {
            auto key = kh_key(seq_locator, k);
            auto val = kh_val(seq_locator, k);

            if (val->size() > 1)
            {
                std::set<uint64_t> genome_ids;
                for(auto id : (*val)){
                    genome_ids.insert(id);
                }
                kh_value(seq_locator, k) = new std::vector<uint64_t>;
                for(auto id : genome_ids) {
                    kh_value(seq_locator, k)->push_back(id);
                    (*index)[id].second++;
                }
            }
        }
    }
}

void printSeqHashMap(khash_t(vector_u64)* seq_locator) {
    uint64_t MAX_HASH = 0;
    uint64_t MIN_HASH = UINT64_MAX;

    std::ofstream seq_locator_file("seq.hashmap");
    for (khiter_t k = kh_begin(seq_locator); k != kh_end(seq_locator); ++k)
    {
      if (kh_exist(seq_locator, k))
      {
        auto key = kh_key(seq_locator, k);
        auto val = kh_value(seq_locator, k);
        seq_locator_file << key << " ";
        for (auto mem : *val)
        {
          seq_locator_file << mem << " ";
          if(mem > MAX_HASH) {
              MAX_HASH = mem;
          }
          if(mem < MIN_HASH) {
              MIN_HASH = mem;
          }
        }
        seq_locator_file << "\n";
      }
    }
    seq_locator_file.close();
    std::ofstream min_max_seq("seq.minmax");
    min_max_seq << "Min value: " << MIN_HASH << std::endl;
    min_max_seq << "Max value: " << MAX_HASH << std::endl;
}

void printNewIndex(Index* index) {
    std::ofstream new_index("MinCE_to_NCBI.index");
    for (int i = 0; i < (*index).size(); i++)
    {
        new_index << i << " " << (*index)[i].first << " " << (*index)[i].second << "\n";
    }
    new_index.close();
}

int main(int argc, char** argv)
{
    std::string path_file = argv[argc-1];
    std::vector<std::string> sequence_paths = readPathFile(path_file); // Paths to all .sequences files

    khash_t(vector_u64)* seq_locator = kh_init(vector_u64); // Initialize our seq_locator hash table
    std::map<std::string, uint64_t> index_mapping; // Initialize map from NCBI_name to MinCE_id
    Index index = read_index("indices", &index_mapping); // Read old indices file to vector of pair<string, int>
    for(auto path : sequence_paths) {
        addSeqToHashMap(path, seq_locator, &index_mapping); // Open .sequences files and log { hash -> MinCE_id } to seq_locator
    }
    // Turn instances of {hash -> 1, 2, 1, 3, 1} to {hash -> 1,2,3} and
    // log number of sequences connected to MinCE_id to our Index vector
    removeDuplicatesAndLogNoOfSeq(seq_locator, &index); 
    printSeqHashMap(seq_locator); // Print seq_locator to file
    printNewIndex(&index); // Print new index to file
    return(0);
}