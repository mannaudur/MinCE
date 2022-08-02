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
        "  -D          Dynamic value of N, based on distance within clique N=(5-dist)*D.\n"
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
    uint16_t N = 20;
    uint16_t D = 0;

    std::pair <size_t,size_t> tsv_dims;
    int opt;
    while ((opt = getopt(argc, argv, "d:k:N:D:f:")) != -1)
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
            case 'D':
                D = atoi(optarg);
                break;
            case 'f':
                batch_file = optarg;
                break;
        }
    }

    std::ofstream failed_cliques("failed_cliques.txt");
    failed_cliques << "Clique\t\tMembers\n" << endl;
    failed_cliques.close();

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
                if(D) {
                    std::size_t inner_dist_pos = clique_path.find_first_of("_")+1;
                    std::string inner_dist = clique_path.substr(inner_dist_pos, 1);
                    N = max((5-stoi(inner_dist))*D,5);
                }
                tsv_dims = make_bit_matrix_from_clique(k, clique_path);
                runSequenceFind(k, N, tsv_dims.first, tsv_dims.second, clique_path);
                std::cout << clique_path << " done" << std::endl;
            }
            fs.close();
        }
}