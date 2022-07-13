#include <bifrost/ColoredCDBG.hpp>
#include <bifrost/CompactedDBG.hpp>

std::pair <size_t,size_t> make_bit_matrix_from_clique(size_t k, std::string clique_file) {
    
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
    

    ColoredCDBG<> dbg(k);
    dbg.buildGraph(opt);
    dbg.buildColors(opt);

    std::ofstream outfile(clique_file + ".fasta");
    size_t i = 0;
    for (const auto& um : dbg) {
        outfile << ">" << i++ << std::endl;
        outfile << um.referenceUnitigToString() << std::endl;
    }
    outfile.close();  

    opt.ratio_kmers = 0.8; // !!! Should this be higher?

    std::vector<std::string> query_files;
    query_files.push_back(clique_file + ".fasta");
    dbg.search(query_files, clique_file, opt.ratio_kmers, false, opt.nb_threads, opt.verbose);
    std::pair <size_t,size_t> tsv_dims;

    tsv_dims = std::make_pair(opt.filename_ref_in.size(),dbg.size());
    return(tsv_dims);
}

std::pair <size_t,size_t> make_bit_matrix_from_fastX(
    size_t k, 
    std::vector<std::string> fastX_files, 
    std::string fastX_name,
    std::vector<std::string> clique_files,
    bool batch) 
    {
    CCDBG_Build_opt opt;
    opt.nb_threads = 4;
    opt.verbose = false;
    opt.outputColors = true;
    
    opt.filename_ref_in = fastX_files;
    
    ColoredCDBG<> dbg(k);
    dbg.buildGraph(opt);
    dbg.buildColors(opt);

    opt.ratio_kmers = 1.0; // !!! Should this be higher?

    
    dbg.search(clique_files, fastX_name, opt.ratio_kmers, false, opt.nb_threads, opt.verbose);
    std::pair <size_t,size_t> tsv_dims;

    tsv_dims = std::make_pair(opt.filename_ref_in.size(),dbg.size());
    return(tsv_dims);
}