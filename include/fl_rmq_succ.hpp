#pragma once

#include <limits>
#include <cstdlib>
#include <climits>
#include <iostream>
#include <cassert>
#include <numeric>
#include <stack>
#include <tuple>
#include <cmath>

#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/find_l2_flat_with.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#include "bp_utils.hpp"
#include "fl_rmq.hpp"

template<typename K, typename Range, typename Pos,
 typename Floating = float, size_t Epsilon = 256>
class ExcessFLRMQ : public FLRMQ<K, Range, Pos, Floating, Epsilon> {

public:

    ExcessFLRMQ() = default;

    explicit ExcessFLRMQ(const std::vector<K> &excess_array) : FLRMQ<K, Range, Pos, Floating, Epsilon>(excess_array) {}

    std::tuple<size_t, size_t, size_t, size_t> query(const size_t i, const size_t j) const {

        const auto k = this->msb(j - i + 1);

        const auto e1 = this->encode(i, i + (1 << k) - 1);
        const auto e2 = this->encode(j - (1 << k) + 1, j);

        const auto first = std::next(this->segments.begin(), this->first_segment[k]);
        const auto last = std::next(this->segments.begin(), this->first_segment[k+1]);

        auto it_s1 = this->segment_for_range(e1, first, last+1);
        auto it_s2 = this->segment_for_range(e2, it_s1, last+1);

        auto [lo1, hi1] = this->reduce_interval(it_s1, e1, this->deltas[k], i, i + (1 << k) - 1);
        auto [lo2, hi2] = this->reduce_interval(it_s2, e2, this->deltas[k], j - (1 << k) + 1, j);

        return std::make_tuple(lo1, hi1, lo2, hi2);
    }

};

template<typename K, typename Range, typename Pos,
 typename Floating = float, size_t Epsilon = 256, size_t Samples = 2048>
class SuccinctFLRMQ {
    
    pasta::BitVector bp;
    pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::BINARY_SEARCH> bp_rs;

    sdsl::int_vector<> exc_samples;
    sdsl::int_vector<> min_exc_samples;
    sdsl::int_vector<> pos_min_exc_samples;

    ExcessFLRMQ<int64_t, Range, Pos, Floating, Epsilon> excess_array_rmq;

    bool reversed;
    size_t n;

public:

    SuccinctFLRMQ() = default;

    explicit SuccinctFLRMQ(const std::vector<K> &data) : n(data.size()) {

        const size_t bp_size = 2 * n + 2;

        bp = pasta::BitVector(bp_size, 0);

        if(build_cartesian_tree<true>(data) < build_cartesian_tree<false>(data)) {
            build_cartesian_tree<true>(data, true);
            reversed = true;
        } else {
            build_cartesian_tree<false>(data, true);
            reversed = false;
        }

        // compute the excess array
        // todo: build it while computing the cartesian tree

        std::vector<int64_t> excess_array;
        excess_array.reserve(bp_size);
        int64_t excess = 0;

        for(auto i = 0; i < bp_size; ++i) {
            excess += bp[i] ? 1 : -1;
            excess_array.push_back(excess);
        }

        sample_excess_array(excess_array);

        excess_array_rmq = ExcessFLRMQ<int64_t, Range, Pos, Floating, Epsilon>(excess_array);

        bp_rs = pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES, 
                                        pasta::FindL2FlatWith::BINARY_SEARCH>(bp);
    }

    [[nodiscard]] size_t query(const size_t i, const size_t j) const {
        size_t m_i = (reversed) ? mirror(i) : i;
        size_t m_j = (reversed) ? mirror(j) : j;

        if(reversed) std::swap(m_i, m_j);

        const size_t bp_i = bp_rs.select1(m_i + 2) - 1;
        const size_t bp_j = bp_rs.select1(m_j + 2);

        const auto [l1, h1, l2, h2] = excess_array_rmq.query(bp_i, bp_j);

        int64_t exc_l1 = 2 * bp_rs.rank1(l1) - l1;
        int64_t exc_l2 = 2 * bp_rs.rank1(l2) - l2;

        // rightmost +-1 rmq using the sampled excess array and precomputed tables

        const auto [pos_min_exc1, min_exc1] = min_excess(exc_l1, l1, h1);
        const auto [pos_min_exc2, min_exc2] = min_excess(exc_l2, l2, h2);

        // notice: the rank operation can be avoided
        // but it requires a lot of operations, don't know
        // if it make sense in practice.
        auto array_pos = ((min_exc1 < min_exc2) ? bp_rs.rank1(pos_min_exc1 + 1) : bp_rs.rank1(pos_min_exc2 + 1)) - 1;

        return (reversed) ? mirror(array_pos) : array_pos;
    }

    inline size_t size() const {
        return excess_array_rmq.size() 
                + (bp_rs.space_usage() * CHAR_BIT) 
                + exc_samples.bit_size()
                + min_exc_samples.bit_size()
                + pos_min_exc_samples.bit_size()
                + bp.size() + sizeof(size_t) + sizeof(bool);
    }

private:

    inline void sample_excess_array(const std::vector<int64_t> &excess_array) {
        const size_t excess_size = excess_array.size();
        const size_t samples = (excess_size / Samples) + 2;

        exc_samples = sdsl::int_vector<>(samples, 0);
        min_exc_samples = sdsl::int_vector<>(samples, 0);
        pos_min_exc_samples = sdsl::int_vector<>(samples, 0);

        for(size_t i = 0; i < excess_size; i += Samples) {
            int64_t min_exc = excess_array[i];
            int64_t pos_min_exc = 0;

            for(size_t j = 1; i + j < excess_size && j < Samples; ++j) {
                if(excess_array[i + j] <= min_exc) {
                    min_exc = excess_array[i + j];
                    pos_min_exc = j;
                }
            }

            exc_samples[i / Samples] = (i + Samples - 1 < excess_size) ? excess_array[i + Samples - 1] : excess_array.back();
            min_exc_samples[i / Samples] = min_exc;
            pos_min_exc_samples[i / Samples] = pos_min_exc;
        }

        sdsl::util::bit_compress(exc_samples);
        sdsl::util::bit_compress(min_exc_samples);
        sdsl::util::bit_compress(pos_min_exc_samples);
    }

    inline size_t mirror(const size_t i) const {
        return n - i - 1;
    }

    inline std::pair<size_t, int64_t> min_excess(const int64_t exc_l, const size_t l, const size_t h) const {
        
        if(h - l < Samples) [[unlikely]] {
            return min_excess_scan(bp, exc_l, l, h);
        }

        const size_t start = l / Samples;
        const size_t end = h / Samples;

        int64_t pos_min_exc = pos_min_exc_samples[start], min_exc = min_exc_samples[start];
        size_t sample_idx = start;

        for(size_t k = start + 1; k <= end; ++k) {
            if(min_exc_samples[k] <= min_exc) {
                min_exc = min_exc_samples[k];
                pos_min_exc = pos_min_exc_samples[k];
                sample_idx = k;
            }
        }

        pos_min_exc = pos_min_exc + (sample_idx * Samples);

        if(pos_min_exc >= l && pos_min_exc <= h) [[likely]] {
            return std::make_pair(pos_min_exc, min_exc);
        } else [[unlikely]] {
            
            // notice (start + 1) * Samples - 1 is never >= bp.size() because of the initial check
            std::tie(pos_min_exc, min_exc) = min_excess_scan(bp, exc_l, l, (start + 1) * Samples - 1);

            // can be avoided

            for(size_t k = start + 1; k < end; ++k) {
                if(min_exc_samples[k] <= min_exc) {
                    min_exc = min_exc_samples[k];
                    pos_min_exc = pos_min_exc_samples[k] + (k * Samples);
                    sample_idx = k;
                }
            }

            const auto [last_min_pos, last_min_exc] = min_excess_scan(bp, exc_samples[end - 1], end * Samples, h);
            pos_min_exc = (last_min_exc <= min_exc) ? last_min_pos : pos_min_exc;
            min_exc = (last_min_exc <= min_exc) ? last_min_exc : min_exc;

            return std::make_pair(pos_min_exc, min_exc);
        }
    }

    template<bool reverse>
    inline bool compare(const K a, const K b) const {
        if constexpr (reverse)
            return a <= b;
        else 
            return a < b;
    }

    template<bool reverse>
    int64_t build_cartesian_tree(const std::vector<K> &data, bool write = false) {
        int64_t max_excess = 0, curr_excess = 0;

        if(data.size() > 0) [[likely]] {
            int64_t curr = (reverse) ? data.size() - 1 : 0;
            size_t bp_curr = 0;
            std::stack<K> s;
            
            s.push(std::numeric_limits<K>::min());

            if(write) {
                bp[bp_curr] = 1;
                bp_curr++;
            }

            while((!reverse && curr < data.size()) || (reverse && curr >= 0)) {

                while(compare<reverse>(data[curr], s.top()) && s.size() > 1) {
                    s.pop();
                    bp_curr++;
                    curr_excess--;
                }

                curr_excess++;

                if(write) {
                    bp[bp_curr] = 1;
                    bp_curr++;
                }

                if(curr_excess > max_excess) max_excess = curr_excess;

                s.push(data[curr]);

                curr += (reverse) ? -1 : 1;
            }
        }

        return max_excess;
    }

};