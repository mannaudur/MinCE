#include "MinceTools.hpp"
#include <getopt.h>

#include <chrono>

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.fasta> | <in.fastq>) ...\n\n"
        "Options:\n"
        
        "  -d <path>   Destination directory for sketch(s).\n"
        "  -t          Threshold t/5000 to display as result [default: 4950].\n"
        "  -c          Candidate set limit [default: 3].\n"
        "  -n          Name of output files.\n"
        "  -b          Run batch file (text file with paths to multiple fastq files, 1 per line).\n";
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
    uint16_t k = 31;
    size_t c = 3;
    uint16_t t = 4950;
    bool b = false;
    uint64_t max_hash = 6725168124769237367;
    std::string hashmapdir = "";
    std::vector<std::string> files;

    int opt;
    while ((opt = getopt(argc, argv, "d:n:c:t:b:")) != -1)
    {
        switch (opt)
        {
            case 'd':
                hashmapdir = optarg;
                if (hashmapdir.back() != '/')
                    hashmapdir += '/';
                break;
            case 'b':
                b = true;
                break;
            case 't':
                t = atoi(optarg);
                break;
            case 'c':
                c = atoi(optarg);
                break;
        }
    }
    char* filename = argv[argc-1];
    std::cout << "\n" << "Mincing " << filename << "...\n" << endl;
    Kmer::set_k(k);

    if(b) {
        std::fstream fs(filename, std::ios::in);
        std::string subfile;
        while (std::getline(fs, subfile)) {
            files.push_back(subfile);
        }
        fs.close();
    } else {
        files.push_back(filename);
    }


    std::vector<Result> results = mince_fastq(filename, files, hashmapdir, max_hash, k, c);
    printResultsToConsole(t, filename, results);
    printResultsToFile(t, filename, results);
    writeTSV(filename, results);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Total runtime: " << duration.count() << " seconds\n" << std::endl;
    return 0;
}