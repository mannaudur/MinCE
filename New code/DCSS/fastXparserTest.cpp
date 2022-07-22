#include <bifrost/CompactedDBG.hpp>

int main(int argc, char* argv[]) {
    std::string fasta_file = argv[argc-1];
    std::vector<std::string> vec {fasta_file};
    char read;
    size_t read_len;
    string seq;
    size_t seq_len;
    unsigned int file_id;
    char* qual = NULL;
    //FastqFile f(vec);
    FileParser f(vec);
    int i = 0;
    std::cout << "hi" << i++ << std::endl;
    std::string DNA_sequence;
    size_t k = 31;
    Kmer sequence;

    while(f.read(seq,seq_len)) {
        std::cout << i++ << std::endl;
        std::cout << seq.substr(0,31) << std::endl;
    }
    
    /* while (f.read_next(&read, &read_len, seq, &seq_len, &file_id, qual) != -1) {
        std::cout << std::stoull(&read) << std::endl;
        std::cout << seq.substr(0,31) << std::endl;   
    } */
    return(0);
}