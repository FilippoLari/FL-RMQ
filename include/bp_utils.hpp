#pragma once

#include <cstdlib>
#include <cstdint>
#include <limits>
#include <utility>

#include <pasta/bit_vector/bit_vector.hpp>

/**
 * Static precomputed table that gives the cumulative excess
 * of any block of 8 bits.
 */
constexpr int8_t exc_byte[] = {
    -8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,
    -6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
    -6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
    -4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
    -6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
    -4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
    -4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
    -2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
    -6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,
    -4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
    -4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
    -2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
    -4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,
    -2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
    -2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,
    0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8,
};

/**
 * Static precomputed table that gives the minimum excess
 * of any block of 8 bits.
 */
constexpr int8_t min_exc_byte[] = {
    -8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,
    -6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,1,
    -6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,1,
    -4,-2,-2,0,-2,0,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,1,
    -4,-2,-2,0,-2,0,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -5,-3,-3,-1,-3,-1,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -4,-2,-2,0,-2,0,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -7,-5,-5,-3,-5,-3,-3,-1,-5,-3,-3,-1,-3,-1,-1,1,
    -5,-3,-3,-1,-3,-1,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -5,-3,-3,-1,-3,-1,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -4,-2,-2,0,-2,0,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,1,
    -4,-2,-2,0,-2,0,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -5,-3,-3,-1,-3,-1,-1,1,-3,-1,-1,1,-2,0,-1,1,
    -4,-2,-2,0,-2,0,-1,1,-3,-1,-1,1,-2,0,-1,1,
};

/**
 * Static precomputed table that gives the position of the rightmost minimum excess
 * of any block of 8 bits.
 */
constexpr int8_t pos_min_exc_byte[] = {
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,0,0,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,0,0,
    7,7,7,7,7,7,0,0,2,2,2,2,1,1,0,0,
    7,7,7,7,7,7,7,7,7,7,7,7,7,7,0,0,
    7,7,7,7,7,7,0,0,2,2,2,2,1,1,0,0,
    4,4,4,4,4,4,4,4,4,4,4,4,1,1,0,0,
    3,3,3,3,3,3,0,0,2,2,2,2,1,1,0,0,
    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
    6,6,6,6,6,6,6,6,6,6,6,6,1,1,0,0,
    6,6,6,6,6,6,6,6,6,6,6,6,1,1,0,0,
    3,3,3,3,3,3,0,0,2,2,2,2,1,1,0,0,
    5,5,5,5,5,5,5,5,5,5,5,5,5,5,0,0,
    5,5,5,5,5,5,0,0,2,2,2,2,1,1,0,0,
    4,4,4,4,4,4,4,4,4,4,4,4,1,1,0,0,
    3,3,3,3,3,3,0,0,2,2,2,2,1,1,0,0,
};

/**
 * Computes the rightmost minimum excess and its position inside a 64-bit word.
 * Notice how exc, min_exc, and min_exc_idx are modified.
 * 
 * Credits: https://github.com/ot/succinct/blob/master/bp_vector.cpp
 * 
 * @param word the 64-bit word
 * @param exc the excess before word
 * @param word_start the starting index of the word inside the bit vector
 * @param min_exc the current minimum excess
 * @param min_exc_idx the current minimum excess position
 */
inline void min_excess_word(const uint64_t word, int64_t &exc, const size_t word_start,
                                                     int64_t &min_exc, size_t &min_exc_idx) {
    int64_t min_byte_exc = std::numeric_limits<int64_t>::max();
    size_t min_byte_idx = 0;

    for(size_t i = 0; i < 8; ++i) {
        const size_t shift = i * 8;
        const size_t byte = (word >> shift) & 0xFF;

        const int64_t curr_min = exc + min_exc_byte[byte];

        min_byte_exc = (curr_min <= min_byte_exc) ? curr_min : min_byte_exc;
        min_byte_idx = (curr_min <= min_byte_exc) ? i : min_byte_idx;

        exc += exc_byte[byte];
    }

    if(min_byte_exc <= min_exc) {
        min_exc = min_byte_exc;
        const size_t shift = min_byte_idx * 8;
        min_exc_idx = word_start + shift + pos_min_exc_byte[(word >> shift) & 0xFF];
    }
}

/**
 * Computes the rightmost minimum excess and its position inside a given range [i,j].
 * 
 * Credits: https://github.com/ot/succinct/blob/master/bp_vector.cpp
 * 
 * @param bv the bit vector
 * @param exc_i the excess of [0,i)
 * @param i the left extreme of the range
 * @param j the right extreme of the range
 * @return a pair consisting of the position of the rightmost minimum excess
 *          and its value 
 */
inline std::pair<size_t, int64_t> min_excess(const pasta::BitVector &bv, int64_t exc_i,
                                                 const size_t i, const size_t j) {
    if (i == j) [[unlikely]] {
        const int64_t excess = bv[i] ? 1 : -1;
        return std::make_pair(i, excess);
    }

    const size_t w_i = i / 64;
    const size_t w_j = j / 64;

    int64_t min_exc = std::numeric_limits<int64_t>::max();
    int64_t curr_exc = exc_i;
    size_t min_exc_idx = i;

    const size_t shift_len_i = (i % 64);
    const size_t range_len = std::min<size_t>(j - i + 1, 64 - shift_len_i);
    const size_t shifted_w_i = bv.data(w_i) >> shift_len_i;
    const size_t padded_w_i = (range_len == 64) ? shifted_w_i : shifted_w_i | (~0ULL << range_len);

    min_excess_word(padded_w_i, curr_exc, i, min_exc, min_exc_idx);

    // todo: [[likely]] or [[unlikely]], check how often does it happen
    if (w_i == w_j) {
        return std::make_pair(min_exc_idx, min_exc);
    }

    curr_exc -= shift_len_i;

    for(size_t w = w_i + 1; w < w_j; ++w) {
        min_excess_word(bv.data(w), curr_exc, w * 64, min_exc, min_exc_idx);
    }

    const size_t shift_len_j = (j % 64) + 1; // includes j
    const size_t padded_w_j = (shift_len_j == 64) ? bv.data(w_j) : bv.data(w_j) | (~0ULL << shift_len_j); 

    min_excess_word(padded_w_j, curr_exc, w_j * 64, min_exc, min_exc_idx);

    return std::make_pair(min_exc_idx, min_exc);
}