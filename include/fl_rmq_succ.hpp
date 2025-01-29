#pragma once

#include <unordered_map>
#include <iostream>
#include <cstdlib>
#include <climits>
#include <cassert>
#include <numeric>
#include <limits>
#include <stack>
#include <tuple>
#include <cmath>

#include <pasta/bit_vector/support/find_l2_flat_with.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#include <pasta/bit_vector/bit_vector.hpp>

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#include "sparse_table.hpp"
#include "bp_utils.hpp"
#include "fl_rmq.hpp"

template<typename K, typename Range, typename Pos,
 typename Floating = float, size_t Epsilon = 4096>
class ExcessFLRMQ : public FLRMQ<K, Range, Pos, Floating, 0, Epsilon> {

public:

    ExcessFLRMQ() = default;

    explicit ExcessFLRMQ(const std::vector<K> &excess_array) : FLRMQ<K, Range, Pos, Floating, 0, Epsilon>(excess_array) {}

    std::tuple<size_t, size_t, size_t, size_t> query(const size_t i, const size_t j) const {

        const auto k = sdsl::bits::hi(j - i + 1);

        const auto e1 = this->encode(i, i + (1 << k) - 1);
        const auto e2 = this->encode(j - (1 << k) + 1, j);

        const auto first = std::next(this->ranges.begin(), this->first_segment[k]);
        const auto last = std::next(this->ranges.begin(), this->first_segment[k+1]);

        auto it_s1 = this->segment_for_range(e1, first, last+1);
        auto it_s2 = this->segment_for_range(e2, it_s1, last+1);

        auto [lo1, hi1] = this->reduce_interval(std::distance(this->ranges.begin(), it_s1),
                                                 e1, this->deltas[k], i, i + (1 << k) - 1);
        auto [lo2, hi2] = this->reduce_interval(std::distance(this->ranges.begin(), it_s2),
                                                 e2, this->deltas[k], j - (1 << k) + 1, j);

        return std::make_tuple(lo1, hi1, lo2, hi2);
    }

};

template<typename K, typename Range, typename Pos,
 typename Floating = float, size_t Epsilon = 4096,
 size_t ExcSamples = 2048, size_t DataSamples = 0>
class SuccinctFLRMQ {
    static_assert(ExcSamples > 0);
    static_assert(ExcSamples <= Epsilon);

    pasta::BitVector bp;
    pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::BINARY_SEARCH> bp_rs;

    sdsl::int_vector<> exc_samples;
    sdsl::int_vector<> min_exc_samples;
    sdsl::int_vector<> pos_min_exc_samples;

    ExcessFLRMQ<int64_t, Range, Pos, Floating, Epsilon> excess_array_rmq;

    sdsl::int_vector<> min_data_samples;
    sdsl::int_vector<> pos_min_data_samples;

    SparseTable<sdsl::int_vector<>> top_level_st;

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

        if constexpr (DataSamples > 0) {
            sample_data(data);
            top_level_st = SparseTable<sdsl::int_vector<>>(&min_data_samples);
        }
    }

    [[nodiscard]] size_t query(const size_t i, const size_t j) const {
        if constexpr (DataSamples > 0) {
            const size_t ds_i = i / DataSamples;
            const size_t ds_j = j / DataSamples;

            const size_t min_sample = top_level_st.query(ds_i, ds_j);
            const size_t min_idx = pos_min_data_samples[min_sample] + min_sample * DataSamples;

            if(i <= min_idx && min_idx <= j) return min_idx;
        }
        
        size_t m_i = (reversed) ? mirror(i) : i;
        size_t m_j = (reversed) ? mirror(j) : j;

        if(reversed) std::swap(m_i, m_j);

        const size_t bp_i = bp_rs.select1(m_i + 2) - 1;
        const size_t bp_j = bp_rs.select1(m_j + 2);

        // linear scan with precomputed tables and sampled
        // excess array if the range is small

        if(bp_j - bp_i + 1 <= 2 * Epsilon + 2) {
            const auto [pos_min, _] = min_excess(bp_i, bp_j);
            const auto array_pos = bp_rs.rank1(pos_min + 1) - 1;
            return (reversed) ? mirror(array_pos) : array_pos;
        }

        // rightmost +-1 rmq on the reduced intervals using
        // the precomputed tables and sampled excess array 

        const auto [l1, h1, l2, h2] = excess_array_rmq.query(bp_i, bp_j);

        const auto [pos_min_exc1, min_exc1] = min_excess(l1, h1);
        const auto [pos_min_exc2, min_exc2] = min_excess(l2, h2);

        // notice: the rank operation can be avoided
        // but it requires a lot of operations, don't know
        // if it make sense in practice.
        auto array_pos = ((min_exc1 < min_exc2) ? bp_rs.rank1(pos_min_exc1 + 1) : bp_rs.rank1(pos_min_exc2 + 1)) - 1;

        return (reversed) ? mirror(array_pos) : array_pos;
    }

    inline size_t size() const {
        const auto data_samples_size = (DataSamples > 0) ? min_data_samples.bit_size() 
                                                            + pos_min_data_samples.bit_size() 
                                                            + top_level_st.size() : 0;
        return excess_array_rmq.size()
                + data_samples_size
                + precomputed_tables_size()
                + (bp_rs.space_usage() * CHAR_BIT) 
                + exc_samples.bit_size()
                + min_exc_samples.bit_size()
                + pos_min_exc_samples.bit_size()
                + bp.size() + (sizeof(size_t) + sizeof(bool)) * CHAR_BIT;
    }

private:

    void sample_data(const std::vector<K> &data) {
        auto samples = n / DataSamples + ((n % DataSamples) != 0);

        std::vector<K> tmp_min_data_samples(samples);
        min_data_samples = sdsl::int_vector<>(samples);
        pos_min_data_samples = sdsl::int_vector<>(samples);

        for(auto i = 0; i < samples; ++i) {
            auto min_idx = i * DataSamples;
            for(auto j = min_idx; j < (i+1) * DataSamples && j < n; ++j) {
                if(data[j] < data[min_idx]) min_idx = j; // notice: leftmost minimum here
            }
            tmp_min_data_samples[i] = data[min_idx];
            pos_min_data_samples[i] = min_idx - i * DataSamples;
        }

        // since data can contain arbitrarily large values
        // we remap its content into [1, n / DataSamples].
        // Notice that this does not change the result of any query

        std::vector<K> sorted_samples = tmp_min_data_samples;
        std::sort(sorted_samples.begin(), sorted_samples.end());

        std::unordered_map<K, size_t> rank_map;
        size_t rank = 0;

        for(const auto &sample : sorted_samples) {
            if(rank_map.find(sample) == rank_map.end())
                rank_map[sample] = rank++;
        }

        for(auto i = 0; i < tmp_min_data_samples.size(); ++i)
            min_data_samples[i] = rank_map[tmp_min_data_samples[i]];

        sdsl::util::bit_compress(min_data_samples);
        sdsl::util::bit_compress(pos_min_data_samples);
    }

    void sample_excess_array(const std::vector<int64_t> &excess_array) {
        const size_t excess_size = excess_array.size();
        const size_t samples = (excess_size / ExcSamples) + 2;

        exc_samples = sdsl::int_vector<>(samples, 0);
        min_exc_samples = sdsl::int_vector<>(samples, 0);
        pos_min_exc_samples = sdsl::int_vector<>(samples, 0);

        for(size_t i = 0; i < excess_size; i += ExcSamples) {
            int64_t min_exc = excess_array[i];
            int64_t pos_min_exc = 0;

            for(size_t j = 1; i + j < excess_size && j < ExcSamples; ++j) {
                if(excess_array[i + j] <= min_exc) {
                    min_exc = excess_array[i + j];
                    pos_min_exc = j;
                }
            }

            exc_samples[i / ExcSamples] = (i + ExcSamples - 1 < excess_size) ? excess_array[i + ExcSamples - 1] : excess_array.back();
            min_exc_samples[i / ExcSamples] = min_exc;
            pos_min_exc_samples[i / ExcSamples] = pos_min_exc;
        }

        sdsl::util::bit_compress(exc_samples);
        sdsl::util::bit_compress(min_exc_samples);
        sdsl::util::bit_compress(pos_min_exc_samples);
    }

    inline size_t mirror(const size_t i) const {
        return n - i - 1;
    }

    inline std::pair<size_t, int64_t> min_excess(const size_t l, const size_t h) const {
        
        if(h - l + 1 < ExcSamples) [[unlikely]] {
            const auto exc_l = 2 * bp_rs.rank1(l) - l;
            return min_excess_scan(bp, exc_l, l, h);
        }

        const size_t start = l / ExcSamples;
        const size_t end = h / ExcSamples;

        int64_t pos_min_exc = pos_min_exc_samples[start], min_exc = min_exc_samples[start];
        size_t sample_idx = start;

        for(size_t k = start + 1; k <= end; ++k) {
            if(min_exc_samples[k] <= min_exc) {
                min_exc = min_exc_samples[k];
                pos_min_exc = pos_min_exc_samples[k];
                sample_idx = k;
            }
        }

        pos_min_exc = pos_min_exc + (sample_idx * ExcSamples);

        if(pos_min_exc >= l && pos_min_exc <= h) [[likely]] {
            return std::make_pair(pos_min_exc, min_exc);
        } else [[unlikely]] {
            
            const auto exc_l = 2 * bp_rs.rank1(l) - l;

            // notice (start + 1) * ExcSamples - 1 is never >= bp.size() because of the initial check
            std::tie(pos_min_exc, min_exc) = min_excess_scan(bp, exc_l, l, (start + 1) * ExcSamples - 1);

            // todo: can be avoided
            for(size_t k = start + 1; k < end; ++k) {
                if(min_exc_samples[k] <= min_exc) {
                    min_exc = min_exc_samples[k];
                    pos_min_exc = pos_min_exc_samples[k] + (k * ExcSamples);
                    sample_idx = k;
                }
            }

            const auto [last_min_pos, last_min_exc] = min_excess_scan(bp, exc_samples[end - 1], end * ExcSamples, h);
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