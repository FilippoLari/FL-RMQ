#pragma once

#include <cstdlib>
#include <cstdint>
#include <limits>
#include <utility>

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
 * Given a 64-bit word return the byte containing the minimum
 * excess and the cumulative excess value of the word.
 */
inline std::pair<size_t, int64_t> min_excess_word(uint64_t w) {

    int64_t curr_min = std::numeric_limits<int64_t>::max();
    int64_t curr_excess = 0;
    size_t min_byte_idx = 0;

    for(size_t i = 0; i < 8; ++i) {
        size_t shift = i * 8;
        size_t byte = (w >> shift) & 0xFF;

        int64_t min_exc = curr_excess + min_exc_byte[byte];

        curr_min = (min_exc <= curr_min) ? min_exc : curr_min;
        min_byte_idx = (min_exc <= curr_min) ? i : min_byte_idx;

        curr_excess += exc_byte[byte];
    }

    return std::make_pair(min_byte_idx, curr_excess);
}