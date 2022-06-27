#include "khash.h"
#include "Sketch.hpp"
#include "UnionFind.hpp"
#include <getopt.h>
#include <cstdint>
#include <algorithm>
#define BYTE_FILE 0
#define CLUSTER_LIMIT_ERROR 5

using Pair = std::pair<uint64_t, uint64_t>;
const bool cmp(const Pair x, const Pair y)
{
    return x.second > y.second;
}

KHASH_MAP_INIT_INT64(vector_u64, std::vector<uint64_t>*);
KHASH_MAP_INIT_INT64(u64, uint64_t);

khash_t(vector_u64)* make_hash_locator(std::vector<Sketch>& sketch_list)
{
    int ret;
    khiter_t k;

    khash_t(vector_u64)* hash_locator = kh_init(vector_u64);

    for (uint64_t i = 0; i < sketch_list.size(); i++)
    {
        Sketch& sketch = sketch_list[i];
        for (uint64_t hash : sketch.min_hash)
        {
            k = kh_get(vector_u64, hash_locator, hash);

            if (k == kh_end(hash_locator))
            {
                k = kh_put(vector_u64, hash_locator, hash, &ret);
                kh_value(hash_locator, k) = new std::vector<uint64_t>;
            }

            kh_value(hash_locator, k)->push_back(i);
        }
    }

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

void make_clusters(
    UnionFind& uf,
    const std::vector<Sketch>& sketch_list,
    khash_t(vector_u64)* hash_locator,
    const uint64_t limit,
    bool names)
{
    int ret;
    khiter_t k;
    std::vector<std::vector<uint64_t>> connections;

    // UnionFind uf(sketch_list.size());

    for (uint64_t i = 0; i < sketch_list.size(); i++)
    {
        // Indices of sketches and number of mutual hash values.
        khash_t(u64)* mutual = kh_init(u64);

        for (auto hash : sketch_list[i].min_hash)
        {
            // Indices of sketches where hash appears.
            k = kh_get(vector_u64, hash_locator, hash);
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
        std::vector<uint64_t> connects_to;
        for (k = kh_begin(mutual); k != kh_end(mutual); ++k)
        {
            if (kh_exist(mutual, k))
            {
                const auto j = kh_key(mutual, k);
                const auto c = kh_value(mutual, k);

                if (c >= limit)
                {
                    connects_to.push_back(j);
                    if (uf.find(i) != uf.find(j)) 
                    {
                        uf.merge(i, j);
                    }
                }
            }
        }
        connections.push_back(connects_to);

        kh_destroy(u64, mutual);
    }
    std::ofstream connect_log("connect_log.txt");
    for (uint64_t i = 0; i < sketch_list.size(); i++)
    {
        if (names) {
            connect_log << sketch_list[i].fastx_filename;
            for (int j = 0; j < connections[i].size(); j++) 
            {
                if (sketch_list[connections[i][j]].fastx_filename != sketch_list[i].fastx_filename) {
                    connect_log << " " << sketch_list[connections[i][j]].fastx_filename;
                }
            }
        } else {
            connect_log << i;
            for (int j = 0; j < connections[i].size(); j++) 
            {
                if(connections[i][j] != i) {
                    connect_log << " " << connections[i][j];
                }
            }
        }
        connect_log << "\n";
    }
    connect_log.close();

    khash_t(vector_u64)* clusters = kh_init(vector_u64);

    for (int x = 0; x < uf.size(); x++)
    {
        const int parent = uf.find(x);

        k = kh_get(vector_u64, clusters, parent);

        if (k == kh_end(clusters))
        {
            k = kh_put(vector_u64, clusters, parent, &ret);
            kh_value(clusters, k) = new std::vector<uint64_t>;
        }

        kh_value(clusters, k)->push_back(x);
    }
}

void usage()
{
    static char const s[] = "Usage: make_graph [options] <file>\n\n"
        "Options:\n"
        "   -l <u64>    Mininum of mutual k-mers [default: 990/1000].\n"
        "   -r          Rep sketch path.\n"
        "   -i          Info file name.\n"
        "   -n          Print file names instead of index.\n"
        "   -h          Show this screen.\n";
    std::printf("%s\n", s);
}

int main(int argc, char** argv)
{
    uint64_t limit = 990;
    std::string rep_path = "";
    std::string info_file = "";
    bool names = false;

    int option;
    while ((option = getopt(argc, argv, "l:r:i:n:h:")) != -1)
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
            case 'n':
                names = true;
                break;
            case 'h':
                usage();
                exit(0);
        }
    }

    std::vector<Sketch> sketch_list;
    std::vector<std::string> fnames = read(argv[optind]);
    sketch_list.reserve(fnames.size());
    for (auto fname : fnames)
    {
        sketch_list.push_back(Sketch::read(fname.c_str()));
    }

    auto hash_locator = make_hash_locator(sketch_list);

    UnionFind uf(sketch_list.size());
    make_clusters(uf, sketch_list, hash_locator, limit, names);
    return 0;
}
