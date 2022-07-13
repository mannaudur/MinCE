#include <bifrost/ColoredCDBG.hpp>
#include <bifrost/CompactedDBG.hpp>
#include <fstream>

int main(int argc, char* argv[]) {

    CCDBG_Build_opt opt;
    opt.nb_threads = 4;
    opt.verbose = false;
    opt.outputColors = true;

    std::string file_prefix = argv[1];
    std::ifstream infile(file_prefix + ".clique");
    std::string filename;
    while (infile >> filename) {
        opt.filename_ref_in.push_back(filename);
    }
    infile.close();
    

    ColoredCDBG<> dbg(31);
    dbg.buildGraph(opt);
    dbg.buildColors(opt);

    std::ofstream outfile(file_prefix + ".fasta");
    size_t i = 0;
    for (const auto& um : dbg) {
        outfile << ">" << i++ << std::endl;
        outfile << um.referenceUnitigToString() << std::endl;
    }
    outfile.close();  

    std::vector<std::string> query_files;
    query_files.push_back(file_prefix + ".fasta");
    dbg.search(query_files, file_prefix, 0.8, false, opt.nb_threads, opt.verbose);
}