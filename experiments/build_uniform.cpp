#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <random>

#include "block_decomposition_hb.hpp"
#include "hybrid_rmq.hpp"
#include "fl_rmq.hpp"

std::vector<int> generate(int n, int lower_bound, int upper_bound) {
    // Check if the range is valid
    if (upper_bound - lower_bound + 1 < n) {
        throw std::invalid_argument("Range is too small to generate unique elements.");
    }

    std::unordered_set<int> unique_numbers;
    std::mt19937 gen(3);
    std::uniform_int_distribution<> dis(lower_bound, upper_bound);

    while (unique_numbers.size() < n) {
        int number = dis(gen);
        unique_numbers.insert(number);
    }

    // Convert set to vector
    std::vector<int> random_vector(unique_numbers.begin(), unique_numbers.end());
    return random_vector;
}

int main(void) {
    
    constexpr size_t eps = 256;
    constexpr size_t threshold = 10000;

    //const std::vector<size_t> sizes = {10000, 100000, 1000000, 10000000};
    const std::vector<size_t> sizes = {1000000};

    std::cout << "eps: " << eps << std::endl;

    for(const auto &size : sizes) {
        std::cout << "Testing " << size << " keys" << std::endl;
        
        std::vector<int> data = generate(size, 1, 1000000000);

        FLRMQ<int, int64_t, int64_t, float, eps> fl_rmq(data);
        block_decomposition_hb<int, int64_t> block_dec(data, static_cast<size_t>(ceil(pow(size, 0.33))));
        hybrid_rmq<int, int64_t, int64_t, float, threshold, eps> hb_rmq(data);

        std::cout << "number of segments: " << fl_rmq.segment_count() << std::endl;
        
        std::cout << "bits per element : " << double(fl_rmq.size()) / double(size) << std::endl;
        std::cout << "bits per element block decomp. : " << double(block_dec.size()) / double(size) << std::endl;
        std::cout << "bits per element hybrid : " << double(hb_rmq.size()) / double(size) << std::endl;

        std::cout << "==========" << std::endl;
    }

    return 0;
}