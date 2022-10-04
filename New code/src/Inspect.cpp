#include <vector>
#include <string>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <algorithm>

std::vector<std::string> splitString(std::string in, char sep) {
    std::vector<std::string> r;
    for (auto p = in.begin();; ++p) {
        auto q = p;
        p = std::find(p, in.end(), sep);
        r.emplace_back(q, p);
        if (p == in.end())
            return r;
    }
}

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] <results.tsv> ...\n\n"
        "Parameters x, t, r, l work in conjunction.\n"
        "Example: '-x 10 -l 3' will output top 10 results having no more than 3 sequences missing.\n\n"
        "Options:\n"
        
        "  -p <path>   Path to fasta files. Easier just to edit the local path permanently in 'Inspect.cpp'.\n"
        "  -n          Name of result-file [default: infile+'.paths'].\n"
        "  -x          Output top-x of results.\n"
        "  -t          Threshold t/5000 to output result.\n"
        "  -r <float>  Ratio of found sequences to output result.\n"
        "  -l          Minimum number of sequences not found to output result.\n"
        "  -c          Print results to console.\n";
    printf(s, name);
}


int main(int argc, char** argv)
{
    if (argc == 1)
    {
        print_usage(argv[0]);
        exit(1);
    }

    std::string path_to_fastas = "/home/thorhallur/newMinCE207/fasta/";
    std::string suffix = "_genomic.fna.gz";
    std::string outname = "";
    size_t x = 0;
    size_t t = 0;
    float r = 0.0;
    int l = 100;
    bool c = false;

    int opt;
    while ((opt = getopt(argc, argv, "p:n:x:t:r:l:c:")) != -1)
    {
        switch (opt)
        {
            case 'd':
                path_to_fastas = optarg;
                if (path_to_fastas.back() != '/')
                    path_to_fastas += '/';
                break;
            case 'n':
                outname = optarg;
                break;
            case 'x':
                x = atoi(optarg);
                break;
            case 't':
                t = atoi(optarg);
                break;
            case 'r':
                r = atof(optarg);
                break;
            case 'l':
                l = atoi(optarg);
                break;
            case 'c':
                c = true;
                break;
        }
    }
    if(!(x || t)) {
        std::cout << "A parameter value for either -x or -t is required." << std::endl;
        return(0);
    }
    std::string filename = argv[argc-1];
    if(outname == "") {
        outname = filename + ".paths";
    }

    std::ifstream tsv_data(filename);
    std::string line;
    std::string name;
    size_t sketch;
    int seq;
    int total_seq;
    float ratio;
    std::vector<std::string> result_line;

    std::ofstream outfile(outname);
    if(c) std::cout << std::endl;
    while(std::getline(tsv_data,line) && (x > 0)) {
        result_line = splitString(line, '\t');
        name = result_line[0];
        sketch = stoi(result_line[1]);
        seq = stoi(result_line[2]);
        total_seq = stoi(result_line[3]);
        ratio = (seq/total_seq);
        if(total_seq < 0) {
            ratio = 1.0;
        }
        if(sketch >= t && ratio >= r && (total_seq - seq) <= l) {
            outfile << path_to_fastas + name + suffix << std::endl;
            if(c) std::cout << path_to_fastas + name + suffix << " ";
            x--;
        }
        if(sketch < t) {
            break;
        }
    }
    if(c) std::cout << std::endl << std::endl;
    outfile.close(); 
    return(0);
}