#ifndef TRIANGULATE_HPP
#define TRIANGULATE_HPP
#include "khash.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>

#include <chrono>


KHASH_MAP_INIT_INT64(vec, std::vector<uint64_t>*);
KHASH_MAP_INIT_INT64(u64, uint64_t);

using HashLocator = khash_t(vec)*;
using Indices = std::vector<std::pair<std::string, int> >;
using Pair = std::pair<uint64_t, uint64_t>;
using Mutual = std::vector<Pair>;

HashLocator read_hash_locator(const char *hash_locator_filename)
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

Indices read_indices(const char *indices_filename)
{
    int index;
    std::string genome;
    std::string atom;
    std::fstream fs(indices_filename, std::ios::in);
    std::string line;
    Indices indices;
    while (std::getline(fs, line))
    {
        std::stringstream ss(line);
        ss >> index >> genome >> atom;
        if (atom == "NULL")
            indices.push_back(std::make_pair(genome, -1));
        else
            indices.push_back(std::make_pair(genome, atoi(atom.c_str())));
    }
    return indices;
}

Mutual find_mutual(uint64_t *min_hash, size_t n)
{
    HashLocator hash_locator = read_hash_locator("../../hash_locator.nosync");
    int ret;
    khint64_t k;
    khint64_t h;
    khash_t(u64) *mutual = kh_init(u64);
    for (int i = 0; i < n; i++)
    {
        k = kh_get(vec, hash_locator, min_hash[i]);
        if (k != kh_end(hash_locator))
        {
            auto indx = kh_value(hash_locator, k);
            for (auto j : *indx)
            {
                h = kh_get(u64, mutual, j);
                if (h != kh_end(mutual))
                {
                    kh_value(mutual, h) += 1;
                }
                else
                {
                    h = kh_put(u64, mutual, j, &ret);
                    kh_value(mutual, h) = 1;
                }
            }
        }
    }
    Mutual mut;
    for (k = kh_begin(mutual); k != kh_end(mutual); ++k)
    {
        if (kh_exist(mutual, k))
            mut.push_back(std::make_pair(kh_key(mutual, k), kh_value(mutual, k)));
    }
    kh_destroy(u64, mutual);

    auto cmp = [](const Pair& x, const Pair& y) -> const bool {
        return x.second > y.second;
    };
    std::sort(mut.begin(), mut.end(), cmp);

    return mut;
}

struct Results
{
    uint64_t id;
    uint64_t mutual;
    int atom;
    char genome[1024];
};

std::vector<Results> get_results(uint64_t *min_hash, size_t min_hash_size, uint16_t threshold)
{
    std::vector<Results> res;

    Mutual mutual = find_mutual(min_hash, min_hash_size);
    Indices indices = read_indices("../../indices.nosync");

    for (size_t i = 0; i < mutual.size(); ++i)
    {
        Results res_mut;
        if(mutual[i].second < threshold)
            break;
        res_mut.id = mutual[i].first;
        res_mut.mutual = mutual[i].second;
        memset(res_mut.genome, 0, sizeof(res_mut.genome));
        strcpy(res_mut.genome, indices[res_mut.id].first.c_str());
        res_mut.atom = indices[res_mut.id].second;
        res.push_back(res_mut);
    }
    return res;
}
#endif
