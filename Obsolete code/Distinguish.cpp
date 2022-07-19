#include <bifrost/ColoredCDBG.hpp>
#include <bifrost/CompactedDBG.hpp>
#include "BitMatrix.hpp"
#include "khash.h"
#include <iomanip>

KHASH_MAP_INIT_INT(u16t, std::vector<uint16_t>*);
using Results = std::pair<std::pair<string,string>, std::pair<int,int>>;

void print_usage(const char *name)
{
    static char const s[] =
        "\nUsage: %s [options] (<in.txt> | <in.fastq>) ...\n\n"
        "Options:\n"
        
        "  -T          Threshold found/available to consider as result.\n"
        "  -t          How many kmers member can miss and still qualify for result (works in conjunction with -T.\n"
        "  -n          Name of results file.\n"
        "  -k          Size of sequences. Will default to closest higher match of 31+(n*32) [default: 31]\n"
        "  -c          .sequences file to read (when -s) or .txt file with paths to multiple .sequences files.\n"
        "  -s          Run single .sequences file\n"
        "  -f          Run a batch of fastq files as input (.txt file with paths).\n";
    printf(s, name);
}

std::vector<Results> readTsvToResults(
    std::string tsv_path, 
    khash_t(u16t)* res_log)
{
    int ret;
    khiter_t k;

    std::ifstream tsv_data(tsv_path+".tsv");
    std::string line;
    std::string member;
    std::string clique;
    std::map<std::string,int> IDsToInt;
    std::map<std::string,std::string> IDsToClique;
    int keep_track = 0;
    int memberID;
    int found;

    std::getline(tsv_data,line); // Don't want first line

    while(std::getline(tsv_data,line))
    {
        std::size_t IDEnd = line.find("_genomic");
        std::size_t CliqueStart = line.find_first_of(';');
        std::size_t CliqueEnd = line.find_last_of(';');

        member = line.substr(1, IDEnd-1);
        clique = line.substr(CliqueStart+1,CliqueEnd-CliqueStart-1);
        found = stoi(line.substr(line.size()-1));
        
        if(IDsToInt.count(member) == 0) {
            IDsToInt[member] = keep_track;
            IDsToClique[member] = clique;
            keep_track++;
        }

        memberID = IDsToInt[member];
        k = kh_get(u16t, res_log, memberID); // Load hits for each member in hash table

        if (k == kh_end(res_log))
        {
            k = kh_put(u16t, res_log, memberID, &ret);
            kh_value(res_log, k) = new std::vector<uint16_t>;
            kh_value(res_log, k)->push_back(0);
            kh_value(res_log, k)->push_back(0);
        }
        auto temp_res = (*kh_value(res_log, k));
        temp_res[1]++;
        temp_res[0]+=found;
        (*kh_value(res_log, k)) = temp_res;
    }

    std::vector<Results> results;
    map<string, int>::iterator itr;
    for (itr = IDsToInt.begin(); itr != IDsToInt.end(); ++itr) {
        k = kh_get(u16t, res_log, (*itr).second);
        if (k != kh_end(res_log)) {
            auto val = (*kh_value(res_log,k));
            auto origin = std::make_pair((*itr).first,IDsToClique[(*itr).first]);
            auto statistics = std::make_pair(val[0],val[1]);
            results.push_back(std::make_pair(origin, statistics));
        }
    }
    std::sort(results.begin(), results.end(),
          [](const Results& a, const Results& b) {
        return a.second.first > b.second.first;
    });
    return(results);
}

void printResults(int T, int t, std::vector<Results> results, std::vector<std::pair<std::string, std::string>> incomplete_results) {
    int found;
    int total;
    std::string NCBI_id;
    std::string clique_id;
    float ratio;
    std::vector<Results> final_results;

    std::cout << "\n NCBI identifier of genome: \t\t|\tClique ID\t| \t  Sequence search results:\t  |\n";
    std::cout << "----------------------------------------|-----------------------|-----------------------------------------|" << std::endl;
    for(auto res : results) {
        std::string space1 = "\t\t";
        NCBI_id = res.first.first;
        if(NCBI_id.size() > 30) {
            space1 = "\t";
        }
        clique_id = res.first.second;
        std::string space2 = "\t\t";
        if(clique_id.size() > 7) {
            space2 = "\t";
        }
        found = res.second.first;
        total = res.second.second;
        ratio = float(found)/float(total);
        
        if(ratio > T || total - found < t+1) {
            std::cout << "                                 \t|                 \t|                                  \t  |" << std::endl;
            std::cout << ' ' << NCBI_id << space1 << "|\t" << clique_id << space2 << "|\t     " << found << '/' << total << "  \t ";
            std::cout << std::fixed << std::setprecision(1) << (ratio)*100 << "%\t\t  |" << std::endl;
        }
    }
    for(auto inc_res : incomplete_results) {
        std::string space1 = "\t\t";
        NCBI_id = inc_res.first;
        if(NCBI_id.size() > 30) {
            space1 = "\t";
        }
        clique_id = inc_res.second;
        std::string space2 = "\t\t";
        if(clique_id.size() > 7) {
            space2 = "\t";
        }

        std::cout << "                                 \t|                 \t|                                  \t  |" << std::endl;
        std::cout << ' ' << inc_res.first << space1 << "|\t" << inc_res.second << space2 << "|\tNo distinguishing sequences  \t  |" << std::endl;
    }
    std::cout << std::endl << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc == 1)
    {
        print_usage(argv[0]);
        exit(1);
    }

    std::string fastX_file = argv[argc-1];
    uint16_t k = 31;
    float T = 0.5;
    int t = 2;

    std::vector<std::string> clique_files;
    std::vector<std::string> fastX_files;
    std::vector<std::string> incomplete_cliques;
    std::vector<std::pair<std::string, std::string>> incomplete_results;
    std::string clique_file_path;
    std::string fastX_name = "";
    std::pair <size_t,size_t> tsv_dims;
    bool batch = false;
    bool single = false;
    int opt;
    while ((opt = getopt(argc, argv, "n:k:c:s:f:T:")) != -1)
    {
        switch (opt)
        {
            case 'n':
                fastX_name = optarg;
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'c':
                clique_file_path = optarg;
                break;
            case 's':
                single = true;
                break;
            case 'f':
                batch = true;
                break;
            case 'T':
                T = atof(optarg);
                break;
        }
    }
    if(batch) {
        std::fstream fs(fastX_file, std::ios::in);
        std::string filename;
        while (std::getline(fs, filename)) {
            fastX_files.push_back(filename);
        }
    } else {
        fastX_files.push_back(fastX_file);
    }
    if(fastX_name == "") {
        fastX_name = fastX_file;
    }
    if(single) {
        std::size_t PosOfFileExt = clique_file_path.find(".sequences");
        std::string noKmers = clique_file_path.substr(PosOfFileExt-2,2);
        if(noKmers == "NA"){
            incomplete_cliques.push_back(clique_file_path);
        } else {
            clique_files.push_back(clique_file_path);
        }
    } else {
        std::fstream fs(clique_file_path, std::ios::in);
        std::string filename;
        while (std::getline(fs, filename)) {
            std::size_t PosOfFileExt = filename.find(".sequences");
            std::string noKmers = filename.substr(PosOfFileExt-2,2);
            if(noKmers == "NA"){
                incomplete_cliques.push_back(filename);
            } else {
                clique_files.push_back(filename);
            }
        }
    }

    for(auto clique_file : incomplete_cliques) {
        std::fstream fs(clique_file, std::ios::in);
        std::string line;
        while (std::getline(fs, line)) {
            if(line[0] == '>') {
                std::size_t IDEnd = line.find("_genomic");

                std::size_t CliqueEnd = clique_file.find_first_of('.');

                std::string member = line.substr(1, IDEnd-1);
                std::string clique = clique_file.substr(0,CliqueEnd);
                auto res = std::make_pair(member, clique);
                incomplete_results.push_back(res);
            }
        }
    }

    std::cout << "\n ---> Mincing " << fastX_file << " ...\n\n";

    tsv_dims = make_bit_matrix_from_fastX(k, fastX_files, fastX_name, clique_files, batch);
    khash_t(u16t)* res_log = kh_init(u16t);

    auto results = readTsvToResults(fastX_name, res_log);
    printResults(T, t, results, incomplete_results);
}

