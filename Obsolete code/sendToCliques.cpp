#include <bifrost/ColoredCDBG.hpp>
#include <bifrost/CompactedDBG.hpp>
#include "BitMatrix.hpp"
#include "khash.h"
#include <iomanip>

KHASH_MAP_INIT_INT(u16t, std::vector<uint16_t>*);
using Origin = std::pair<string,string>;
using Statistics = std::pair<std::pair<int,int>,std::pair<int,int>>;
using Result = std::pair<Origin, Statistics>;
using Final_results = std::vector<Result>;

std::vector<Result> readTsvToResults(
    std::string tsv_path, 
    Final_results results)
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

    for(size_t i = 0; i < results.size(); i++) 
    {
            IDsToInt[results[i].first.first] = i;
    }

    while(std::getline(tsv_data,line))
    {
        std::size_t IDEnd = line.find("_genomic");
        std::size_t CliqueStart = line.find_first_of(';');
        std::size_t CliqueEnd = line.find_last_of(';');

        member = line.substr(1, IDEnd-1);
        found = stoi(line.substr(line.size()-1));
        
        memberID = IDsToInt[member];

        results[memberID].second.second.first+=found;
        results[memberID].second.second.second++;
    }
    
    return(results);
}

void printResults(int T, int t, Final_results results, std::vector<std::pair<std::string, std::string>> incomplete_results) {
    int found;
    int total;
    int sketch_found;
    std::string NCBI_id;
    std::string clique_id;
    float ratio;

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
        found = res.second.second.first;
        total = res.second.second.second;
        ratio = float(found)/float(total);
        sketch_found = res.second.first.first;
        
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

Final_results distinguish(int T, int t, int k, std::vector<std::string> fastX_files, std::string fastX_name,  Final_results results) {
    std::vector<std::string> clique_files;
    std::string header;
    std::string seq_file;
    for(Result res : results) {
        if(res.first.second != "NULL") {
            seq_file = res.first.second+".sequences";
            std::fstream fs(seq_file, std::ios::in);
            std::getline(fs, header);
            if(header[0] == '<') {
                res.second.first = std::make_pair(0,-1);
            } else {
                clique_files.push_back(res.first.second+".sequences");
            }
        }
    }
    make_bit_matrix_from_fastX(k, fastX_files, fastX_name, clique_files, true);
    results = readTsvToResults(fastX_name, results);
}

