#pragma once

#include <cstdlib>
#include <climits>
#include <numeric>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

template<typename K, typename Pos, size_t BlockSize>
class BlockDecomposition {
    static_assert(std::is_integral_v<K>);
    static_assert(std::is_integral_v<Pos>);
    static_assert(BlockSize > 0);

    sdsl::int_vector<> samples;

    size_t n;

public:

    BlockDecomposition() = default;

    explicit BlockDecomposition(const std::vector<K> &data) : n(data.size()) {
        if(n == 0) [[unlikely]]
            return;

        const size_t blocks = (n / BlockSize) + 2;

        samples = sdsl::int_vector<>(blocks, 0);

        for(size_t i = 0; i < n; i += BlockSize) {

            K min = data[i];
            size_t idx = 0;

            for(size_t j = 1; i + j < n && j < BlockSize; ++j) {
                if(data[i + j] < min) {
                    min = data[i + j];
                    idx = j;
                }
            }

            samples[i / BlockSize] = idx;
        }

        sdsl::util::bit_compress(samples);
    }

    /**
     * Returns the position of the minimum inside the interval [i,j]
     * 
     * @param i the left extreme of the interval
     * @param j the right extreme of the interval
     * @return the position of the minimum inside [i,j] 
     */
    Pos query(const std::vector<K> &data, const size_t i, const size_t j) const {
        if(i==j) [[unlikely]]
            return i;
        
        const auto block_start =  i / BlockSize;
        const auto block_end = j / BlockSize;

        if(j - i + 1 <= 2*BlockSize) [[unlikely]] {
            return linear_scan_min(data, i, j).second;
        }

        size_t first_block_idx = samples[block_start] + block_start * BlockSize;
        K first_min = data[first_block_idx];

        size_t covered_block_idx = samples[block_start + 1] + (block_start + 1) * BlockSize;
        K covered_min = data[covered_block_idx];
        size_t block_idx = covered_block_idx;

        for(size_t k = block_start + 2; k < block_end; ++k) {
            block_idx = samples[k] + k * BlockSize;
            if(data[block_idx] < covered_min) {
                covered_min = data[block_idx];
                covered_block_idx = block_idx;
            }
        }

        size_t last_block_idx = samples[block_end] + block_end * BlockSize;
        K last_min = data[last_block_idx];

        K min = (covered_min < first_min) ? covered_min : first_min;
        size_t idx = (covered_min < first_min) ? covered_block_idx : first_block_idx;
        idx = (last_min < min) ? last_block_idx : idx;

        if(idx >= i && idx <= j) [[likely]] {
            return idx;
        } else [[unlikely]] {
            std::tie(first_min, first_block_idx) = linear_scan_min(data, i, (block_start + 1) * BlockSize - 1);

            idx = (covered_min < first_min) ? covered_block_idx : first_block_idx;
            min = (covered_min < first_min) ? covered_min : first_min;

            std::tie(last_min, last_block_idx) = linear_scan_min(data, block_end * BlockSize, j);

            return (last_min < min) ? last_block_idx : idx;
        }
    }

    inline size_t size() const {
        return samples.bit_size() + (sizeof(n) * CHAR_BIT);
    }

private:

    inline std::pair<K, size_t> linear_scan_min(const std::vector<K> &data, const size_t lo, const size_t hi) const {
        K min = data[lo];
        size_t idx = lo;
        
        for(auto i = lo; i <= hi; ++i)
            if(data[i] < min) {
                min = data[i];
                idx = i;
            }
        
        return std::make_pair(min, idx);
    }
};