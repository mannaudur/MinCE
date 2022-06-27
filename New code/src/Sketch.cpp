#include "Sketch.hpp"
#include <getopt.h>

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.fasta> | <in.fastq>) ...\n\n"
        "Options:\n"
        "  -k    Size of kmers [default: 31].\n"
        "  -c    Candidate set limit [default: 1].\n"
        "  -d    Directory to place sketches.\n"
        "  -s    Size of min hash [default: 5000]. Ignored with -x and "
            "-X options.\n"
        "  -x    Include all hashes that have a value lower than "
            "9999999776999205UL.\n"
        "  -f    Read in txt of fasta file paths to sketch.\n"
        "  -X    Include all hashes that have a value lower than X.\n\n";
    printf(s, name);
}

int main(int argc, char** argv)
{
    if (argc == 1)
    {
        print_usage(argv[0]);
        exit(1);
    }

    uint16_t k = 31;
    uint16_t c = 1;
    uint16_t s = 5000;
    bool x = false;
    uint64_t max_hash;
    std::string dirpath = "";
    std::string batch_file = "";

    int opt;
    while ((opt = getopt(argc, argv, "d:f:k:c:s:x:")) != -1)
    {
        switch (opt)
        {
            case 'd':
                dirpath = optarg;
                if (dirpath.back() != '/')
                    dirpath += '/';
                break;
            case 'f':
                batch_file = optarg;
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'c':
                c = atoi(optarg);
                break;
            case 's':
                s = atoi(optarg);
                break;
            case 'x':
                x = true;
                max_hash = atol(optarg);
                break;
        }
    }

    Kmer::set_k(k);

    if (x)
    {
        if (batch_file == "")
        {
            for (; optind < argc; optind++)
                Sketch::xsketch(argv[optind], max_hash, k, c)
                    .write(dirpath);
        }
        else
        {
            std::fstream fs(batch_file, std::ios::in);
            std::string filename;
            while (std::getline(fs, filename))
            {
                Sketch::xsketch(filename.c_str(), max_hash, k, c)
                    .write(dirpath);
            }
            fs.close();
        }
    }
    else
    {
        if (batch_file == "")
        {
            for (; optind < argc; optind++)
            {
                Sketch _sketch = Sketch::sketch(argv[optind], k, c, s);
                _sketch.write(dirpath);
            }
        }
        else
        {
            std::fstream fs(batch_file, std::ios::in);
            std::string filename;
            while (std::getline(fs, filename))
                Sketch::sketch(filename.c_str()).write(dirpath);
            fs.close();
        }
    }

    return 0;
}
