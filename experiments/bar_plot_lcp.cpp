#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream>
#include <chrono>
#include <filesystem>

#include "fl_rmq.hpp"
#include "utils.hpp"

int main(int argc, char* argv[]) {
    
    if(argc < 2) {
        std::cerr << "Please provide a directory" << std::endl;
        std::cerr << "Usage: ./test_lcp directory" << std::endl;
        return -1;
    }

    std::string ds_path = argv[1];

    // header line for the csv file
    // size is in bytes and time is in seconds 
    // ratio is segments / size and bpe is the number of bits per element
    std::cout << "dataset," << "size," << "eps," << "time," << "space," << "segments," 
              << "ratio," << "bpe" << std::endl;  

    for (const auto& entry : std::filesystem::directory_iterator(ds_path)) {
        const std::string path = entry.path();
        const std::string dataset = entry.path().filename();
        const std::string extension = entry.path().extension();

        if (std::filesystem::is_regular_file(path) && extension != LCP_EXT) {

            std::vector<int64_t> lcp = build_lcp(path);
            const size_t n = lcp.size();

            std::cout << "size: " << n << std::endl;

            auto collect_stats = [&dataset, &lcp, &n] (auto eps) {
            
                auto start = std::chrono::high_resolution_clock::now();

                FLRMQ<int64_t, int64_t, int64_t, float, eps> ds(lcp);

                auto end = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> duration = end - start;

                const double time = duration.count();
                const size_t segments = ds.segment_count();
                const double ratio = double(segments) / double(n);
                const double bpe = double(ds.size()) / double(n); 

                std::cout << dataset << "," << n << "," << eps << "," 
                        << time << "," << ds.size() << "," << segments << "," << ratio << "," 
                        << bpe << std::endl;
            };

            static_for_pw<64, 8192>(collect_stats);
        }
    }

    return 0;
}