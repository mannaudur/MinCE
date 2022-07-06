#include "khash.h"
#include "Sketch.hpp"
#include "UnionFind.hpp"
#include <getopt.h>
#include <cstdint>
#include <algorithm>
#define BYTE_FILE 0
#define CLUSTER_LIMIT_ERROR 5

KHASH_MAP_INIT_INT64(vec, std::vector<uint64_t>*);
KHASH_MAP_INIT_INT64(u64, uint64_t);
using Pair = std::pair<uint64_t, uint64_t>;

const bool cmp(const Pair x, const Pair y)
{
    return x.second > y.second;
}

khash_t(vec)* read_hash_locator(const char *hash_locator_filename)
{
    std::cout << "Loading hash_locator..." << endl;
    int ret;
    khint64_t k;
    uint64_t hash;
    uint64_t index;
    khash_t(vec) *hash_locator = kh_init(vec);
    std::fstream fs(hash_locator_filename, std::ios::in);
    std::string line;
    auto start = std::chrono::high_resolution_clock::now();
    while (std::getline(fs, line))
    {
        std::stringstream ss(line);
        ss >> hash;
        k = kh_put(vec, hash_locator, hash, &ret);
        kh_value(hash_locator, k) = new std::vector<uint64_t>;
        while (ss >> index)
            kh_value(hash_locator, k)->push_back(index);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<chrono::seconds>(stop - start);
    std::cout << "Time for hash_locator setup: " << duration.count() << " seconds\n" << std::endl;
    return hash_locator;
}

std::vector<std::string> read(std::string ifpath)
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

khash_t(vec)* make_clusters(
    UnionFind& uf,
    const std::vector<Sketch>& sketch_list,
    khash_t(vec)* hash_locator,
    const uint64_t limit)
{
    int ret;
    khiter_t k;

    // UnionFind uf(sketch_list.size());

    for (uint64_t i = 0; i < sketch_list.size(); i++)
    {
        // Indices of sketches and number of mutual hash values.
        khash_t(u64)* mutual = kh_init(u64);

        for (auto hash : sketch_list[i].min_hash)
        {
            // Indices of sketches where hash appears.
            k = kh_get(vec, hash_locator, hash);
            std::vector<uint64_t>* sketch_indices = kh_value(hash_locator, k);

            for (auto j : *sketch_indices)
            {
                k = kh_get(u64, mutual, j);

                if (k != kh_end(mutual))
                {
                    kh_value(mutual, k) += 1;
                }
                else
                {
                    k = kh_put(u64, mutual, j, &ret);
                    kh_value(mutual, k) = 1;
                }
            }
        }

        for (k = kh_begin(mutual); k != kh_end(mutual); ++k)
        {
            if (kh_exist(mutual, k))
            {
                const auto j = kh_key(mutual, k);
                const auto c = kh_value(mutual, k);

                if (c > limit && uf.find(i) != uf.find(j))
                {
                    uf.merge(i, j);
                }
            }
        }

        kh_destroy(u64, mutual);
    }

    khash_t(vec)* clusters = kh_init(vec);

    for (int x = 0; x < uf.size(); x++)
    {
        const int parent = uf.find(x);

        k = kh_get(vec, clusters, parent);

        if (k == kh_end(clusters))
        {
            k = kh_put(vec, clusters, parent, &ret);
            kh_value(clusters, k) = new std::vector<uint64_t>;
        }

        kh_value(clusters, k)->push_back(x);
    }

    return clusters;
}

void print_usage(const char *name)
{
    static char const s[] = "Usage: %s [options] <file>\n\n"
        "Options:\n"
        "   -l <u64>    Mininum of mutual k-mers [default: 4980/5000].\n"
        "   -r          hash_locator path [default: in folder]\n"
        "   -i          Info file name\n"
        "   -d          Directory for cliques generated.\n"
        "   -h          Show this screen.\n";
    std::printf("%s\n", s);
}

int main(int argc, char** argv)
{
    if (argc == 1)
    {
        print_usage(argv[0]);
        exit(1);
    }

    uint64_t limit = 4990;
    std::string rep_path = "hash_locator";
    std::string info_file = "clique_log";
    std::string dirpath = "cliques/";

    int option;
    while ((option = getopt(argc, argv, "l:r:i:h:")) != -1)
    {
        switch (option)
        {
            case 'l':
                limit = std::atoi(optarg);
                break;
            case 'r':
                rep_path = optarg;
                break;
            case 'i':
                info_file = optarg;
                break;
            case 'd':
                dirpath = optarg;
                if (dirpath.back() != '/')
                    dirpath += '/';
                break;

        }
    }

    std::vector<Sketch> sketch_list;
    std::vector<std::string> fnames = read(argv[optind]);
    sketch_list.reserve(fnames.size());
    for (auto fname : fnames)
    {
        sketch_list.push_back(Sketch::read(fname.c_str()));
    }
    auto hash_locator = read_hash_locator(rep_path.c_str());
    
    
    UnionFind uf(sketch_list.size());
    auto clusters = make_clusters(uf, sketch_list, hash_locator, limit);
    std::ofstream indices("indices");
    for (int i = 0; i < sketch_list.size(); i++)
    {
        auto parent = uf.find(i);
        khiter_t k = kh_get(vec, clusters, parent);

        auto val = kh_val(clusters, k);
        if (val->size() > 1)
        {
            // Modification required if sketch is member of multiple cliques
            indices << i << " " << sketch_list[i].fastx_filename << " " << parent << "\n";
        }
        else
        {
            indices << i << " " << sketch_list[i].fastx_filename.c_str() << " NULL\n";
        }
    }
    indices.close();
    std::vector<std::vector<uint64_t>> clique_log;
    for (khiter_t k = kh_begin(clusters);
         k != kh_end(clusters);
         ++k)
    {
      if (kh_exist(clusters, k))
      {
        auto key = kh_key(clusters, k);
        auto val = kh_val(clusters, k);

        if (val->size() > 1)
        {
          std::ofstream cliques(dirpath + std::to_string(key) + ".clique");

          for (auto i : *val)
            // i is a reference to the line in indices corresponding to said sketch
            // needed to modify later when cliques become overlapping
            cliques << i << " " << sketch_list[i].fastx_filename << std::endl;

          cliques.close();
          clique_log.push_back( { key, val->size() } );
        }
      }
    }
    // Taken from https://stackoverflow.com/a/14419565
    sort(clique_log.begin(), clique_log.end(),[](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) 
    {
        return a[1] > b[1];
    });

    std::ofstream clique_logger(info_file);
    clique_logger << "Number of cliques: " << clique_log.size() << "\n" << std::endl;
    for (size_t i = 0; i < clique_log.size(); i++) {
        clique_logger << clique_log[i][0] << "\t" << clique_log[i][1] << std::endl;
    }
    clique_logger.close();
}
