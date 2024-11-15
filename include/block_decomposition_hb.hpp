#pragma once

#include <limits>
#include <cstdlib>
#include <climits>
#include <numeric>
#include <vector>

template<typename K, typename Pos>
class block_decomposition_hb {
    static_assert(std::is_integral_v<K>);
    static_assert(std::is_integral_v<Pos>);

    std::vector<Pos> block_minima;
    std::vector<K> data;

    size_t block_size;
    size_t n;

public:

    block_decomposition_hb() = default;

    explicit block_decomposition_hb(const std::vector<K> &data,
     size_t block_size) : data(data), block_size(block_size), n(data.size()) {
        if(n == 0) [[unlikely]]
            return;

        const size_t blocks = (n + block_size - 1) / block_size;

        block_minima.reserve(blocks);

        for(auto i = 0; i < n; i += block_size) {
            K min_value = data[i];
            Pos min_pos = i; 

            for(auto j = 0; (j < block_size && i + j < n); ++j) {
                if(data[i + j] < min_value) {
                    min_value = data[i + j];
                    min_pos = i+j;
                } 
            }

            block_minima.push_back(min_pos);
        }
    }

    /**
     * Returns the position of the minimum inside the interval [i,j]
     * 
     * @param i the left extreme of the interval
     * @param j the right extreme of the interval
     * @return the position of the minimum inside [i,j] 
     */
    Pos query(size_t i, size_t j) const {
        if(i==j) [[unlikely]]
            return i;
        
        const auto first_block =  i / block_size;
        const auto last_block = j / block_size;

        auto min_value = data[i];
        auto min_pos = i;

        if(first_block == last_block) {
            for(auto k = i; k <= j; ++k) {
                if(data[k] < min_value) {
                    min_value = data[k];
                    min_pos = k;
                } 
            }
        } else {
            for(auto k = i; k <= (first_block + 1) * block_size - 1; ++k) {
                if(data[k] < min_value) {
                    min_value = data[k];
                    min_pos = k;
                }
            }

            for(auto k = first_block + 1; k < last_block; ++k) {
                if(data[block_minima[k]] < min_value) {
                    min_value = data[block_minima[k]];
                    min_pos = block_minima[k];
                }
            }

            for(auto k = last_block * block_size; k <= j; ++k) {
                if(data[k] < min_value) {
                    min_value = data[k];
                    min_pos = k;
                }
            }
        }

        return min_pos;
    }

    inline size_t size() const {
        return ((block_minima.size() * sizeof(Pos)) 
                    + sizeof(block_size)
                    + sizeof(n)) * CHAR_BIT;
    }
};