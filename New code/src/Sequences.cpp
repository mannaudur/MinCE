#include "Sequences.hpp"
#include "BitMatrix.hpp"

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.txt> | <in.clique>) ...\n\n"
        "Options:\n"
        
        "  -d <path>   Destination directory for clique files.\n"
        "  -k          Size of sequences. Will default to closest higher match of 31+(n*32) [default: 31]\n"
        "  -N          Number of sequences per genome in clique [default: 10]\n"
        "  -f          Run batch file.\n";
    printf(s, name);
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
    uint16_t k = 31;
    uint16_t N = 10;

    std::pair <size_t,size_t> tsv_dims;
    int opt;
    while ((opt = getopt(argc, argv, "d:k:N:f:")) != -1)
    {
        switch (opt)
        {
            case 'd':
                dirpath = optarg;
                if (dirpath.back() != '/')
                    dirpath += '/';
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'N':
                N = atoi(optarg);
                break;
            case 'f':
                batch_file = optarg;
                break;
        }
    }
    if (batch_file == "")
        {
            tsv_dims = make_bit_matrix_from_clique(k, infile);
            runSequenceFind(k, N, tsv_dims.first, tsv_dims.second, infile);
            std::cout << infile << " done" << std::endl;
        }
        else
        {
            std::fstream fs(infile, std::ios::in);
            std::string clique_path;
            while (std::getline(fs, clique_path))
            {
                tsv_dims = make_bit_matrix_from_clique(k, clique_path);
                runSequenceFind(k, N, tsv_dims.first, tsv_dims.second, clique_path);
                std::cout << clique_path << " done" << std::endl;
            }
            fs.close();
        }
}