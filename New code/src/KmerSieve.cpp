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
        "  -c          .minced file to read (when -s) or .txt file with paths to multiple .minced files.\n"
        "  -s          Run single .minced file\n"
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
        std::size_t PosOfClique = line.find_last_of(";");
        std::size_t PosOfCliqueEnd = line.find_last_of("!");

        member = line.substr(0, PosOfClique-12);
        clique = line.substr(PosOfClique,PosOfCliqueEnd);
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
        clique_files.push_back(clique_file_path);
    } else {
        std::fstream fs(clique_file_path, std::ios::in);
        std::string filename;
        while (std::getline(fs, filename)) {
            clique_files.push_back(filename);
        }
    }

    std::cout << "\nMincing " << fastX_file << " ...\n\n";

    tsv_dims = make_bit_matrix_from_fastX(k, fastX_files, fastX_name, clique_files, batch);
    khash_t(u16t)* res_log = kh_init(u16t);

    auto results = readTsvToResults(fastX_name, res_log);
    int found;
    int total;
    std::string NCBI_id;
    std::string clique_id;
    float ratio;

    std::cout << "\n NCBI identifier of genome: \t |    Clique ID    | \t\t Results:\n";
    std::cout << "---------------------------------|-----------------|---------------------------------" << std::endl;
    std::cout << "                                 |                 |                                 " << std::endl;
    for(auto res : results) {
        NCBI_id = res.first.first;
        clique_id = res.first.second;
        found = res.second.first;
        total = res.second.second;
        ratio = float(found)/float(total);
        if(ratio > T || total - found < t+1) {
            std::cout << ' ' << NCBI_id << "\t | " << clique_id << " |\t     " << found << '/' << total << "  \t ";
            std::cout << std::fixed << std::setprecision(1) << (ratio)*100 << "%" << std::endl;
    std::cout << "                                 |                 |                                 " << std::endl;
        }
    }
}

