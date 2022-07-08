#include <fstream>
#include <getopt.h>
#include <vector>

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.txt> | <in.clique>) ...\n\n"
        "Options:\n"
        
        "  -d <path>   Destination directory for sequence files.\n"
        "  -k          Size of sequences [default: 31]\n"
        "  -N          Number of sequences per genome in clique [default: 10]";
    printf(s, name);
}

int** readBitMatrix(std::string clique_path)
    {
    std::ifstream data(clique_path+".tsv");
    std::string line;
    std::vector<std::vector<std::string>> parsedCsv;
    while(std::getline(data,line))
    {
        // Write my own iterator to convert each line into int
        // The big question is whether I need the dimensions of 
        // the matrix beforehand (pretty likely...)
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<std::string> parsedRow;
        while(std::getline(lineStream,cell,','))
        {
            parsedRow.push_back(cell);
        }

        parsedCsv.push_back(parsedRow);
    }
    }

void chooseSequences(std::string clique_path, uint16_t k, uint16_t N)
{

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
    uint16_t k = 31;
    uint16_t N = 10;

    int opt;
    while ((opt = getopt(argc, argv, "d:k:N:")) != -1)
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
        }
    }


}