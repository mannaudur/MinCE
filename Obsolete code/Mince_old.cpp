#include "Sketch.hpp"
#include "Triangulate.hpp"
#include <getopt.h>

#include <chrono>

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.fasta> | <in.fastq>) ...\n\n"
        "Options:\n"
        
        "  -d <path>   Destination directory for sketch(s).\n"
        "  -x          Megasketch: Include all hashes that have a value lower than "
            "9999999776999205UL.\n"
        "  -t          Threshold t/1000 to display as result [default: 850].\n"
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
    string filename = argv[argc-1];
    std::cout << "\n" << "Mincing " << filename << "...\n" << endl;
    uint16_t k = 31;
    uint16_t c = 3;
    uint16_t s = 1000;
    uint16_t t = 850;
    bool x = false;
    bool p = false;
    uint64_t max_hash;
    std::string dirpath = "";
    //std::string batch_file = "";
    Sketch _sketch;
    std::vector<uint64_t> min_hash;

    // Setup phase before anything really happens
    int opt;
    while ((opt = getopt(argc, argv, "d:f:t:c:x:")) != -1)
    {
        switch (opt)
        {
            case 'd':
                dirpath = optarg;
                if (dirpath.back() != '/')
                    dirpath += '/';
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
            case 'x':
                x = true;
                max_hash = atol(optarg);
                break;
            case 'p':
                p = true;
                break;
        }
    }

    Kmer::set_k(k);

    // Sketch phase, sketch saved to memory
    if (x)
    {
        for (; optind < argc; optind++) 
            {
            _sketch = Sketch::xsketch(argv[optind], max_hash, k, c);
            if(p)
                _sketch.write(dirpath);
            }   
    }
    else
    {
        for (; optind < argc; optind++)
            {
            _sketch = Sketch::sketch(argv[optind], k, c, s);
            if(p)
                _sketch.write(dirpath);
            }
    }
    for (; optind < argc; optind++) {
        std::cout << argv[optind] << std::endl;
    }

    // Load necessary data for distance approximation from new input sketch
        // Might be better to configure the hashlocator as a B-tree instead of a hash-table mapping
        // We should expect some changes in the way distance is calculated and in the way clusters are
        // distributed, f.x. with the advent of overlapping clusters

    // Getting results of close clusters/sketches or multiple hits for a megasketch
    // Size of results saved for future reference
    auto results = get_results(_sketch.min_hash.data(), _sketch.min_hash.size(), t);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<chrono::seconds>(stop - start);
    std::cout << "\n" << "Total duration: " << duration.count() << " seconds\n" << std::endl;

    std::ofstream file;
    file.open(dirpath + get_filename_from_path(filename) + "_mince.txt");
    std::cout << "Results from mincing " << filename << ":\n" << std::endl;
    file << "Results from mincing " << filename << ":\n\n";
    for(size_t i = 0; i < results.size(); i++) {
        std::cout << i+1 << ":\t" << results[i].genome << "\t\t\t" << results[i].mutual << "/" << s << "\t" 
        << "Clique: " << results[i].clique << std::endl;
        file << i+1 << ":\t" << results[i].genome << "\t" << results[i].mutual << "/" << s << "\t" 
        << "Clique: " << results[i].clique << "\n\n";
    }
    file.close();

    return 0;
}