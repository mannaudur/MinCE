#ifndef SKETCH_HPP_
#define SKETCH_HPP_

#include <bifrost/CompactedDBG.hpp>

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <zlib.h>

//KSEQ_INIT(gzFile, gzread);

struct CandidateSet {
    const size_t update(const Kmer);
    void erase(const Kmer);

    private:

    KmerHashTable<uint32_t> candidates;
    KmerHashTable<uint32_t>::iterator it;
};

const size_t CandidateSet::update(const Kmer kmer)
{
    it = candidates.find(kmer);
    return it != candidates.end() ? ++(*it) : candidates.insert(kmer, 1).second;
}

void CandidateSet::erase(const Kmer kmer)
{
    candidates.erase(kmer);
}

std::string get_filename_from_path(std::string path)
{
    size_t last_slash_index = path.find_last_of("\\/");
    if (std::string::npos != last_slash_index)
        path.erase(0, last_slash_index + 1);
    return path;
}

struct Sketch
{
    std::vector<uint64_t> min_hash;
    std::string fastx_filename;
    uint64_t max_hash = 0;
    uint32_t k;
    uint32_t c;
    uint32_t s;

    Sketch() = default;

    void write(std::string dirpath = "")
    {
        std::string extension = max_hash ? ".xsketch" : ".sketch";
        std::ofstream file;
        file.open(dirpath + get_filename_from_path(fastx_filename) + extension);
        file << fastx_filename << '\n';
        file << k << '\n';
        file << c << '\n';
        file << s << '\n';
        for (auto hash : min_hash)
            file << hash << '\n';
        file.close();
    }

    static Sketch read(const char *sketch_filename)
    {
        Sketch _sketch;
        std::fstream fs(sketch_filename, std::ios::in);
        fs >> _sketch.fastx_filename;
        fs >> _sketch.k >> _sketch.c >> _sketch.s;
        std::istream_iterator<uint64_t> start(fs), end;
        _sketch.min_hash.assign(start, end);
        fs.close();
        return _sketch;
    }

    static Sketch sketch(const char *fastx_filename, uint32_t k_ = 31,
            uint32_t c_ = 1, uint32_t s_ = 5000)
    {
        Sketch _sketch;
        _sketch.fastx_filename = fastx_filename;
        _sketch.k = k_;
        _sketch.c = c_;
        _sketch.s = s_;
        Kmer::set_k(k_);

        gzFile fp = gzopen(fastx_filename, "r");
        kseq_t *seq = kseq_init(fp);

        CandidateSet set;
        KmerIterator it_end;
        _sketch.min_hash.reserve(_sketch.s);

        int seq_len;
        while ((seq_len = kseq_read(seq)) >= 0)
        {
            KmerIterator it(seq->seq.s);

            for (; _sketch.min_hash.size() < _sketch.s && it != it_end; ++it)
            {
                const Kmer kmer = it->first.rep();

                if (_sketch.c == set.update(kmer))
                {
                    _sketch.min_hash.push_back(kmer.hash());
                }
            }

            if (_sketch.min_hash.size() == _sketch.s)
                std::make_heap(_sketch.min_hash.begin(), _sketch.min_hash.end());

            for (; it != it_end; ++it)
            {
                const Kmer kmer = it->first.rep();

                if (_sketch.min_hash[0] > kmer.hash())
                {
                    if (_sketch.c == set.update(kmer))
                    {
                        _sketch.min_hash[0] = kmer.hash();
                        sink(_sketch.min_hash);
                    }
                }
            }
        }

        std::sort(_sketch.min_hash.begin(), _sketch.min_hash.end());

        kseq_destroy(seq);
        gzclose(fp);

        return _sketch;
    }

    static Sketch xsketch(const char *fastx_filename,
            uint64_t max_hash, uint32_t k = 31, uint32_t c = 1)
    {
        Sketch sketch;
        sketch.fastx_filename = fastx_filename;
        sketch.k = k;
        sketch.c = c;
        sketch.max_hash = max_hash;
        Kmer::set_k(k);

        gzFile fp = gzopen(fastx_filename, "r");
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

                if (kmer.hash() < max_hash)
                {
                    if (c == set.update(kmer))
                    {
                        sketch.min_hash.push_back(kmer.hash());
                    }
                }
            }
        }
        sketch.s = sketch.min_hash.size();
        std::sort(sketch.min_hash.begin(), sketch.min_hash.end());

        kseq_destroy(seq);
        gzclose(fp);

        return sketch;
    }

    private:

    static void sink(std::vector<uint64_t>& min_hash)
    {
        size_t i = 0;
        size_t j = 1;

        while (j < min_hash.size())
        {
            const size_t temp = j + 1;
            if (temp < min_hash.size() && min_hash[temp] > min_hash[j])
                j = temp;

            if (min_hash[i] > min_hash[j])
                return;

            const uint64_t hash = min_hash[i];
            min_hash[i] = min_hash[j];
            min_hash[j] = hash;

            i = j;
            j = j << 1;
        }
    }
};

#endif
