#include <bifrost/ColoredCDBG.hpp>
#include <bifrost/CompactedDBG.hpp>
#include "run_bifrost.hpp"

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in>) ...\n\n"
        "\nBe advised: In file should be input with extension, i.e. no '.txt' \n\n"
        "Options:\n"
        "  -r    Reference .txt file to parse to bitmatrix.\n"
        "  -d    Directory to place bitmatrix (defult: current directory).\n";
    printf(s, name);
}


int main(int argc, char** argv) 
{
    CCDBG_Build_opt opt;
    opt.build = true;
    int option;
    opt.k = 31;
    std::string dirpath = "";
    while ((option = getopt(argc, argv, "k:d:")) != -1)
    {
        switch (option)
        {
            case 'k':
                opt.k = atoi(optarg);
                break;
            case 'd':
                dirpath = optarg;
                if (dirpath.back() != '/')
                    dirpath += '/';
                break;

        }
    }
    std::string file_prefix = argv[optind];
    opt.filename_ref_in.push_back(file_prefix+".txt");
    opt.prefixFilenameOut = dirpath + file_prefix;
    
    if(!check_ProgramOptions(opt)) {
        std::cout << "Program options check failed" << std::endl;
        return(0);
    }
    //Build colored CDBG
    ColoredCDBG<> ccdbg(opt.k, opt.g); //Initialize ccdbg
    bool success = ccdbg.buildGraph(opt); //Build it
    if (success) success = ccdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose); //...simplify (?)
    if (success) success = ccdbg.buildColors(opt); //Color it
    //if (success) success = ccdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.verbose);

    //Build compacted CDBG
    CompactedDBG<> cdbg(opt.k, opt.g); //Initialize cdbg
    opt.outputGFA = false; //Set output to fasta file
    success = cdbg.build(opt); //Build it
    if (success) success = cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose); //...simplify (?)
    if (success) success = cdbg.write(opt.prefixFilenameOut, opt.nb_threads, opt.outputGFA, opt.verbose); //write to fasta file


    //Query colored GFA with fasta file

    opt.filename_query_in.push_back(opt.prefixFilenameOut + ".fasta"); //Set target query to fasta file from cdbg

    //Query ccdbg with fasta file from cdbg and write results to .tsv file
    success = ccdbg.search(opt.filename_query_in, opt.prefixFilenameOut, opt.ratio_kmers, opt.inexact_search, opt.nb_threads, opt.verbose);
    
    
}



