#include <bifrost/CompactedDBG.hpp>

int main(int argc, char* argv[]) {
    std::string fasta_file = argv[argc-1];
    std::vector<std::string> vec;
    char read;
    size_t read_len;
    string seq;
    size_t seq_len;
    unsigned int file_id;
    char* qual = NULL;
    vec.push_back(fasta_file);
    FastqFile f(vec);
    int i = 0;
    std::cout << "hi" << i++ << std::endl;
    size_t* ret;
    std::string DNA_sequence;
    size_t k = 31;
    Kmer sequence;

    while (f.read_next(&read, &read_len, seq, &seq_len, &file_id, qual)) {
        std::cout << std::stoull(&read) << std::endl;
        std::cout << seq.substr(0,31) << std::endl;
    }

    /*  
    while (f.read_next(s, *ret) ) {
        std::cout << *ret << std::endl;
        std::cout << s.substr(0,31) << std::endl;
    
        DNA_sequence = s.substr(0,min(k, *seq_len));
        seq = Kmer(DNA_sequence.c_str()).rep(); // rep() gets canonical Kmer
        std::cout << atoi(read) << '\n' << seq.toString() << std::endl;

    } */
}