#include <bifrost/ColoredCDBG.hpp>
#include <bifrost/CompactedDBG.hpp>
#include <fstream>

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.txt> | <in.clique>) ...\n\n"
        "Options:\n"
        
        "  -d <path>   Destination directory for sketch(s).\n"
        "  -f          Run batch file.\n";
    printf(s, name);
}

size_t make_bit_matrix(std::string clique_file, std::string dir_path) {
    CCDBG_Build_opt opt;
    opt.nb_threads = 4;
    opt.verbose = false;
    opt.outputColors = true;

    std::ifstream infile(clique_file);
    std::string filename;
    while (infile >> filename) {
        opt.filename_ref_in.push_back(filename);
    }
    infile.close();
    

    ColoredCDBG<> dbg(31);
    dbg.buildGraph(opt);
    dbg.buildColors(opt);

    std::ofstream outfile(clique_file + ".fasta");
    size_t i = 0;
    for (const auto& um : dbg) {
        outfile << ">" << i++ << std::endl;
        outfile << um.referenceUnitigToString() << std::endl;
    }
    outfile.close();  

    std::vector<std::string> query_files;
    query_files.push_back(clique_file + ".fasta");
    dbg.search(query_files, clique_file, 0.8, false, opt.nb_threads, opt.verbose);
    return(dbg.size());
}

int main(int argc, char* argv[]) 
{
    if (argc == 1)
    {
        print_usage(argv[0]);
        exit(1);
    }

    std::string infile = argv[argc-1];
    std::string dirpath = "";
    std::string batch_file = "";

    size_t size_of_sequence_file;
    int opt;
    while ((opt = getopt(argc, argv, "d:f:")) != -1)
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
        }
    }
    if (batch_file == "")
        {
            for (; optind < argc; optind++)
                size_of_sequence_file = make_bit_matrix(infile, dirpath);
                std::cout << size_of_sequence_file << std::endl;
        }
        else
        {
            std::fstream fs(infile, std::ios::in);
            std::string clique_path;
            while (std::getline(fs, clique_path))
            {
                size_of_sequence_file = make_bit_matrix(clique_path, dirpath);
                std::cout << size_of_sequence_file << std::endl;
                std::cout << clique_path << " done" << std::endl;
            }
            fs.close();
        }
}