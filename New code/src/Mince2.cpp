#include "MinceTools.hpp"
#include <getopt.h>

#include <chrono>

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.fasta> | <in.fastq>) ...\n\n"
        "Options:\n"
        
        "  -d <path>   Destination directory for sketch(s).\n"
        "  -a          Run MinCE on a .fasta file (already aligned).\n"
        "  -t          Threshold t/5000 to display as result [default: 4950].\n"
        "  -c          Candidate set limit [default: 3].\n"
        "  -f          Run batch file.\n"
        "  -p          Print sketch(es) to file.\n";
    printf(s, name);
}


int main(int argc, char** argv)
{
    if (argc == 1)
    {
        print_usage(argv[0]);
        exit(1);
    }
    auto start = std::chrono::high_resolution_clock::now();
    char* filename = argv[argc-1];
    std::cout << "\n" << "Mincing " << filename << "...\n" << endl;
    uint16_t k = 31;
    uint16_t c = 3;
    uint16_t s = 5000;
    uint16_t t = 4950;
    bool x = false;
    bool p = false;
    uint64_t max_hash;
    std::string hashmapdir = "";

    // Setup phase before anything really happens
    int opt;
    while ((opt = getopt(argc, argv, "d:f:t:c:x:")) != -1)
    {
        switch (opt)
        {
            case 'd':
                hashmapdir = optarg;
                if (hashmapdir.back() != '/')
                    hashmapdir += '/';
                break;
            /*case 'f':
                batch_file = optarg;
                break;*/
            case 't':
                t = atoi(optarg);
                break;
            case 'c':
                c = atoi(optarg);
                break;
            case 'p':
                p = true;
                break;
        }
    }

    Kmer::set_k(k);
    std::string line;
    std::fstream fs("MAX_HASH.log", std::ios::in);
    std::getline(fs, line);
    max_hash = std::stoull(line);

    std::vector<Result> results = process(filename, hashmapdir, max_hash, k, c, p);
    printResultsToConsole(t, filename, results);
    printResultsToFile(t, filename, filename, results); //BLEGH, FIX THIS

    return 0;
}