#include "Sketch.hpp"
#include <fstream>

size_t shared_kmers(const std::vector<uint64_t> &m1,
                    const std::vector<uint64_t> &m2,
                    int hash_size)
{

    size_t shared = 0;
    size_t i = 0;
    size_t j = 0;

    while (i != hash_size && j != hash_size) {
        if (m1[i] == m2[j]) {
            shared += 1;
            i += 1;
            j += 1;
        } else if (m1[i] < m2[j]) {
            i += 1;
        } else if (m1[i] > m2[j]) {
            j += 1;
        }
    }

    return shared;
}

int main(int argc, char* argv[]) {

    int hash_size = 5000;
    std::string file_prefix = argv[1];
    std::vector<Sketch> sketches;
    std::ifstream infile(file_prefix+".clique");
    std::string filename;
    while (infile >> filename) {
        auto s = Sketch::read(filename.c_str());
        sketches.push_back(s);
        std::cout << s.fastx_filename << " " << s.s << " " << s.min_hash.size() << std::endl;
    }
    infile.close();

    std::ofstream outfile(file_prefix + ".dist");
    for(int i = 0; i < sketches.size(); i++) {
        for(int j = 0; j < sketches.size(); j++) {
            int dist = hash_size - shared_kmers(sketches[i].min_hash, sketches[j].min_hash, hash_size);
            outfile << dist << " ";
        }
        outfile << std::endl;
    }
    outfile << endl;
    for (const auto& sk : sketches) {
        outfile << sk.fastx_filename << std::endl;
    }
    outfile.close();
    return(0);
}