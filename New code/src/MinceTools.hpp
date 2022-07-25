#include <bifrost/CompactedDBG.hpp>
#include "khash.h"
#include "Sketch.hpp"
#include <iomanip>

KHASH_MAP_INIT_INT64(vec, std::vector<uint64_t>*);
KHASH_MAP_INIT_INT(vec32, std::vector<uint32_t>*);

using Index = std::vector<std::pair<std::string, int>>;
using HashLocator = khash_t(vec)*;
using Stats = std::pair<int, std::pair<int, int>>;
using Origin = std::pair<std::string, int>;
using Result = std::pair<Origin, Stats>;


HashLocator read_sketch_hashmap(const char *sketch_hashmap_filename)
{
    std::cout << "Loading sketch.hashmap..." << endl;
    int ret;
    khint64_t k;
    uint64_t hash;
    uint64_t index;
    khash_t(vec) *sketch_hashmap = kh_init(vec);
    std::fstream fs(sketch_hashmap_filename, std::ios::in);
    std::string line;
    auto start = std::chrono::high_resolution_clock::now();
    while (std::getline(fs, line))
    {
        std::stringstream ss(line);
        ss >> hash;
        k = kh_put(vec, sketch_hashmap, hash, &ret);
        kh_value(sketch_hashmap, k) = new std::vector<uint64_t>;
        while (ss >> index)
            kh_value(sketch_hashmap, k)->push_back(index);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time for sketch.hashmap setup: " << duration.count() << " seconds\n" << std::endl;
    return sketch_hashmap;
}

HashLocator read_seq_hashmap(const char *seq_hashmap_filename)
{
    std::cout << "Loading seq.hashmap..." << endl;
    int ret;
    khint64_t k;
    uint64_t hash;
    uint64_t index;
    khash_t(vec) *seq_hashmap = kh_init(vec);
    std::fstream fs(seq_hashmap_filename, std::ios::in);
    std::string line;
    auto start = std::chrono::high_resolution_clock::now();
    while (std::getline(fs, line))
    {
        std::stringstream ss(line);
        ss >> hash;
        k = kh_put(vec, seq_hashmap, hash, &ret);
        kh_value(seq_hashmap, k) = new std::vector<uint64_t>;
        while (ss >> index)
            kh_value(seq_hashmap, k)->push_back(index);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time for seq.hashmap setup: " << duration.count() << " seconds\n" << std::endl;
    return seq_hashmap;
}

Index read_index(const char *index_filename)
{
    std::cout << "Loading index..." << endl;
    int index_id;
    std::string genome;
    std::string NOsequences;
    std::fstream fs(index_filename, std::ios::in);
    std::string line;
    Index index;
    auto start = std::chrono::high_resolution_clock::now();
    while (std::getline(fs, line))
    {
        std::stringstream ss(line);
        ss >> index_id >> genome >> NOsequences;
        index.push_back(std::make_pair(genome, stoi(NOsequences)));
        /// ATH INDEX SHOULD HAVE { ID   -   NCBI_NAME   -    NUMBER_OF_SEQUENCES_CONNECTED } FORMAT
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time for index setup: " << duration.count() << " seconds\n" << std::endl;
    return index;
}

void log_common(khash_t(vec) *hash_locator, khash_t(vec32) *mutual, uint64_t hash, int id) {
    khiter_t k;
    khiter_t h;
    int ret;
    k = kh_get(vec, hash_locator, hash);
    if (k != kh_end(hash_locator))
    {
        auto indx = kh_value(hash_locator, k);
        for (auto j : *indx)
        {
            h = kh_get(vec32, mutual, j);
            if (h != kh_end(mutual))
            {
                (*kh_value(mutual, h))[id] += 1;
            }
            else
            {
                h = kh_put(vec32, mutual, j, &ret);
                kh_value(mutual, h) = new std::vector<uint32_t>;
                kh_value(mutual, h)->push_back(0);
                kh_value(mutual, h)->push_back(0);
                (*kh_value(mutual, h))[id] = 1;
            }
        }
    }
}

std::vector<Result> process(
    const char* filename, 
    std::string hashmapdir, 
    uint64_t max_hash_, 
    uint32_t k_, 
    int c_, 
    bool p)
{
    // Strip sketch from this once I don't need debugging anymore.
    Sketch sketch;
    sketch.fastx_filename = filename;
    sketch.k = k_;
    sketch.c = c_;
    sketch.max_hash = max_hash_;
    Kmer::set_k(k_);

    std::string sketch_hashmap_path = hashmapdir+"sketch.hashmap";
    std::string seq_hashmap_path = hashmapdir+"seq.hashmap";
    HashLocator sketch_hashmap = read_sketch_hashmap(sketch_hashmap_path.c_str());
    HashLocator seq_hashmap = read_seq_hashmap(seq_hashmap_path.c_str());

    auto start = std::chrono::high_resolution_clock::now();
    khash_t(vec32) *mutual = kh_init(vec32);

    gzFile fp = gzopen(filename, "r");
    kseq_t *seq = kseq_init(fp);

    CandidateSet set;
    KmerIterator it_end;

    int seq_len;
    while ((seq_len = kseq_read(seq)) >= 0)
    {
        KmerIterator it(seq->seq.s);

        for (; it != it_end; ++it)
        {
            const Kmer kmer = it->first.rep();
            auto hashmer = kmer.hash();
            if (c_ == set.update(kmer))
            {
                if (hashmer <= max_hash_)
                {
                    log_common(sketch_hashmap, mutual, hashmer, 0);
                    sketch.min_hash.push_back(kmer.hash());
                }
                log_common(seq_hashmap, mutual, hashmer, 1);
            }
        }
    }
    kh_destroy(vec, sketch_hashmap);
    kh_destroy(vec, seq_hashmap);
    kseq_destroy(seq);
    gzclose(fp);
    if(p){
        sketch.s = sketch.min_hash.size();
        std::sort(sketch.min_hash.begin(), sketch.min_hash.end());
        sketch.write();
    }
    std::string index_path = hashmapdir+"MinCE_to_NCBI.index";
    Index index = read_index(index_path.c_str());
    khiter_t k;
    std::vector<Result> results;
    for (k = kh_begin(mutual); k != kh_end(mutual); ++k)
    {
        if (kh_exist(mutual, k)) {
            std::string NCBI_name = index[kh_key(mutual, k)].first; // Log NCBI name of genome
            int NoOfSeq = index[kh_key(mutual, k)].second; // Log number of sequences expected to find for genome

            auto val = kh_value(mutual, k);
            auto sketch_found = (*val)[0]; // Sketch comparison results
            auto seq_found = (*val)[1]; // Sequence comparison results
            
            Origin org = std::make_pair(NCBI_name,kh_key(mutual, k));
            Stats stat = std::make_pair(sketch_found, std::make_pair(seq_found,NoOfSeq));
            
            results.push_back(std::make_pair(org, stat));
        }       
    }
    kh_destroy(vec32, mutual);

    auto cmp1 = [](const Result& x, const Result& y) -> const bool {
        return (x.second.second.second - x.second.second.first) > (y.second.second.second - y.second.second.first);
    };
    std::sort(results.begin(), results.end(), cmp1); // Sorted by number of unique sequences found (or fewest ones missing)

    auto cmp2 = [](const Result& x, const Result& y) -> const bool {
        return x.second.first > y.second.first;
    };
    std::sort(results.begin(), results.end(), cmp2); // Sorted by sketch distance, so we also get members not in cliques
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time spent just mincing: " << duration.count() << " seconds\n" << std::endl;
    return(results);
    }

void printResultsToConsole(int T, std::string infile, std::vector<Result> results) {
    int seq_found;
    int seq_total;
    int sketch_found;
    std::string NCBI_id;
    std::string clique_id;
    float ratio;

    std::cout  << "MinCE results for: " << infile << std::endl;

    std::cout << "\n NCBI identifier of genome: \t\t|   Sketch results:\t| \t  Sequence search results:\t  |\n";
    std::cout << "----------------------------------------|-----------------------|-----------------------------------------|" << std::endl;
    for(auto res : results) {
        sketch_found = res.second.first;
        if(sketch_found < T) {
            break;
        }
        std::string space1 = "\t\t";
        NCBI_id = res.first.first;
        if(NCBI_id.size() > 30) {
            space1 = "\t";
        }
        seq_found = res.second.second.first;
        seq_total = res.second.second.second;
        ratio = 0;
        if(seq_total != 0)
            ratio = float(seq_found)/float(seq_total);
        
        sketch_found = res.second.first;
        
        std::cout << "                                 \t|                 \t|                                  \t  |" << std::endl;
        std::cout << ' ' << NCBI_id << space1 << "|\t" << sketch_found << "/5000\t" << "|\t     " << seq_found << '/' << seq_total << "  \t ";
        std::cout << std::fixed << std::setprecision(1) << (ratio)*100 << "%\t\t  |" << std::endl;
    }
    std::cout << std::endl << std::endl;
}

void printResultsToFile(int T, std::string infile, std::vector<Result> results) {
    int seq_found;
    int seq_total;
    int sketch_found;
    std::string NCBI_id;
    std::string clique_id;
    float ratio;

    std::ofstream outfile(infile + ".minced");

    outfile << infile << std::endl;

    outfile << "\n NCBI identifier of genome: \t\t|   Sketch results:\t| \t  Sequence search results:\t  \n";
    outfile << "------------------------------------|-------------------|---------------------------------" << std::endl;
    for(auto res : results) {
        sketch_found = res.second.first;
        if(sketch_found < T) {
            break;
        }
        std::string space1 = "\t\t";
        NCBI_id = res.first.first;
        if(NCBI_id.size() > 30) {
            space1 = "\t";
        } else if (NCBI_id.size() < 27) {
            space1 = "\t\t\t";
        }
        seq_found = res.second.second.first;
        seq_total = res.second.second.second;
        ratio = float(seq_found)/float(seq_total);
        sketch_found = res.second.first;
        
        outfile << "                                 \t|                 \t|                               " << std::endl;
        outfile << ' ' << NCBI_id << space1 << "|    " << sketch_found << "/5000\t\t" << "|\t     " << seq_found << '/' << seq_total << "  \t\t";
        outfile << std::fixed << std::setprecision(1) << (ratio)*100 << "%" << std::endl;
    }
}