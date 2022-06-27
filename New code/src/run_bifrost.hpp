#include <bifrost/ColoredCDBG.hpp>
#include <bifrost/CompactedDBG.hpp>

bool check_ProgramOptions(CCDBG_Build_opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    auto check_files = [&](vector<string>& v_files) {

        vector<string> files_tmp;

        char* buffer = new char[4096]();

        for (const auto& file : v_files) {

            if (!check_file_exists(file)) {

                cerr << "Error: File " << file << " not found." << endl;
                ret = false;
            }
            else {

                const int format = FileParser::getFileFormat(file.c_str());

                if (format >= 0) files_tmp.push_back(file); // File is FASTA/FASTQ/GFA
                else {

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp != NULL){

                        fclose(fp);

                        ifstream ifs_file_txt(file);
                        istream i_file_txt(ifs_file_txt.rdbuf());

                        size_t i = 0;

                        while (i_file_txt.getline(buffer, 4096).good()){

                            fp = fopen(buffer, "r");

                            if (fp == NULL) {

                                cerr << "Error: Could not open file at line " << i << " in file " << file << " for reading." << endl;
                                ret = false;
                                break;
                            }
                            else {

                                fclose(fp);
                                files_tmp.push_back(string(buffer));
                            }

                            ++i;
                        }

                        if (i_file_txt.fail() && (i == 0)) {

                            cerr << "Error: File " << file << " is neither FASTA, FASTQ nor GFA." << endl;
                            cerr << "If it is a list of files, it is either empty or has a line with >4096 characters." << endl;
                            ret = false;
                        }

                        ifs_file_txt.close();
                    }
                    else {

                        cerr << "Error: Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                }
            }
        }

        v_files = move(files_tmp);

        delete[] buffer;
    };

    // Check general parameters

    if (!opt.build && !opt.update && !opt.query){

        cerr << "Error: No command selected (can be 'build' or 'update' or 'query')." << endl;
        ret = false;
    }

    if (opt.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
        ret = false;
    }

    if (opt.k <= 2){

        cerr << "Error: Length k of k-mers cannot be less than 3." << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed " << (MAX_KMER_SIZE - 1) << "." << endl;
        cerr << "To enable a larger k, recompile Bifrost with the appropriate MAX_KMER_SIZE variable." << endl;
        ret = false;
    }

    if (opt.g == 0){

        cerr << "Error: Length m of minimizers cannot be equal to 0." << endl;
        ret = false;
    }

    if ((opt.g >= 0) && (opt.g > opt.k - 2)) {

        cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
        ret = false;
    }

    if (opt.query){  // Check param. command build

        if (opt.prefixFilenameOut.length() == 0) {

            cerr << "Error: No output filename prefix given." << endl;
            ret = false;
        }
        else {

            const string out = opt.prefixFilenameOut + ".tsv";

            FILE* fp = fopen(out.c_str(), "w");

            if (fp == NULL) {

                cerr << "Error: Could not open file for writing output of query in TSV format: " << out << "." << endl;
                ret = false;
            }
            else {

                fclose(fp);
                if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
            }
        }

        if (opt.filename_query_in.size() == 0) {

            cerr << "Error: Missing input query files." << endl;
            ret = false;
        }
        else check_files(opt.filename_query_in);

        if ((opt.ratio_kmers < 0.0) || (opt.ratio_kmers > 1.0)) {

            cerr << "Error: Ratio of k-mers from queries that must occur in the graph cannot be less than 0.0 or more than 1.0 (" << opt.ratio_kmers << ")." << endl;
            ret = false;
        }

        if (opt.g > opt.k - 2){

            cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
            ret = false;
        }
    }
    else {

        if (opt.prefixFilenameOut.length() == 0) {

            cerr << "Error: No output filename prefix given." << endl;
            ret = false;
        }
        else {

            const string out = opt.prefixFilenameOut + (opt.outputGFA ? ".gfa" : ".fasta");

            FILE* fp = fopen(out.c_str(), "w");

            if (fp == NULL) {

                cerr << "Error: Could not open file for writing output graph in GFA format: " << out << "." << endl;
                ret = false;
            }
            else {

                fclose(fp);
                if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
            }
        }

        if ((opt.filename_seq_in.size() + opt.filename_ref_in.size()) == 0) {

            cerr << "Error: Missing input files." << endl;
            ret = false;
        }
        else {

            check_files(opt.filename_seq_in);
            check_files(opt.filename_ref_in);
        }
    }

    if (opt.build){ // Check param. command build

        if (opt.outFilenameBBF.length() != 0){

            FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

            if (fp == NULL) {

                cerr << "Error: Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
                ret = false;
            }
            else {

                fclose(fp);

                if (remove(opt.outFilenameBBF.c_str()) != 0){

                    cerr << "Error: Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
                }
            }
        }

        if (opt.inFilenameBBF.length() != 0){

            if (check_file_exists(opt.inFilenameBBF)){

                FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "Error: Input Blocked Bloom filter " << opt.inFilenameBBF << " file does not exist." << endl;
                ret = false;
            }
        }
    }

    if (opt.update || opt.query){

        if (opt.filename_graph_in.length() == 0){

            cerr << "Error: No graph file was provided in input." << endl;
            ret = false;
        }
        else if (!check_file_exists(opt.filename_graph_in)){

            cerr << "Error: The graph file does not exist." << endl;
            ret = false;
        }
        else {

            FILE* fp = fopen(opt.filename_graph_in.c_str(), "r");

            if (fp == NULL) {

                cerr << "Error: Could not read input graph file " << opt.filename_graph_in << "." << endl;
                ret = false;
            }
            else fclose(fp);
        }

        if (opt.filename_colors_in.length() != 0){

            if (!check_file_exists(opt.filename_colors_in)){

                cerr << "Error: The input color file does not exist." << endl;
                ret = false;
            }
            else {

                FILE* fp = fopen(opt.filename_colors_in.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input color file " << opt.filename_colors_in << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
        }
    }

    return ret;
}