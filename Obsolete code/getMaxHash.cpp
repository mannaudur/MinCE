#include <fstream>
#include <sstream>
#include <iostream>

int main(int argc, char* argv[])
{
    uint64_t MAX_HASH = 0;
    std::fstream fs("sketch.hashmap", std::ios::in);
    std::string line;
    uint64_t hash;
    while (std::getline(fs, line))
    {
        std::stringstream ss(line);
        ss >> hash;
        if(hash > MAX_HASH)
            MAX_HASH = hash;
    }
    std::cout << MAX_HASH << std::endl;
    return(0);
}