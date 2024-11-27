#pragma once

#include <limits>
#include <cstdlib>
#include <climits>
#include <iostream>
#include <cassert>
#include <numeric>
#include <stack>
#include <tuple>

#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/find_l2_flat_with.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>

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
 typename Floating = float, size_t Epsilon = 256>
class SuccinctFLRMQ {
    
    pasta::BitVector bp;
    pasta::FlatRankSelect<pasta::OptimizedFor::ONE_QUERIES,
                            pasta::FindL2FlatWith::BINARY_SEARCH> bp_rs;

    ExcessFLRMQ<int64_t, Range, Pos, Floating, Epsilon> excess_array_rmq;

    bool reversed;
    size_t n;

public:

    SuccinctFLRMQ() = default;

    explicit SuccinctFLRMQ(const std::vector<K> &data) : n(data.size()) {

        bp = pasta::BitVector(2 * n + 2);

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
        excess_array.reserve(2 * n + 2);
        int64_t excess = 0;

        for(auto i = 0; i < 2 * n + 2; ++i) {
            if(bp[i]) excess++;
            else excess--;
            excess_array.push_back(excess);
        }

        /*int64_t max = 0;
        for(auto &e : excess_array)
            if(e > max) max = e;*/

        //std::cout << "max excess: " << max << std::endl;

        /*for(auto i = 0; i < 2 * n + 2; ++i)
            if(bp[i]) std::cout << "(";
            else std::cout << ")";
        std::cout << std::endl;

        for(auto &e : excess_array)
            std::cout << e << " ";
        std::cout << std::endl;*/

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

        const auto [lo1, hi1, lo2, hi2] = excess_array_rmq.query(bp_i, bp_j);

        // naive linear scans

        int64_t min_excess1 = std::numeric_limits<int64_t>::max();
        int64_t min_excess2 = std::numeric_limits<int64_t>::max(); 
        int64_t excess1_pos = lo1, excess2_pos = lo2;
        int64_t excess1 = 0, excess2 = 0;

        // rightmost +-1 rmq

        for(auto i = lo1; i <= hi1; ++i) {
            if(bp[i]) ++excess1;
            else --excess1;
            min_excess1 = (excess1 <= min_excess1) ? excess1 : min_excess1;
        }

        for(auto i = lo2; i <= hi2; ++i) {
            if(bp[i]) ++excess2;
            else --excess2;
            min_excess2 = (excess2 <= min_excess2) ? excess2 : min_excess2;
        }

        auto array_pos = ((min_excess1 < min_excess2) ? bp_rs.rank1(excess1_pos) : bp_rs.rank1(excess2_pos)) - 1;

        return (reversed) ? mirror(array_pos) : array_pos;
    }

    inline size_t size() const {
        std::cout << "pm RMQ: " << excess_array_rmq.size() << " bp: " << bp.size() << " rs_bp: " << bp_rs.space_usage() * CHAR_BIT << std::endl;
        return excess_array_rmq.size() 
                + (bp_rs.space_usage() * CHAR_BIT) 
                + bp.size() + sizeof(size_t) + sizeof(bool);
    }

private:

    inline size_t mirror(const size_t i) const {
        return n - i - 1;
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

        //std::cout << "reverse: " << reverse << " ";

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

        //std::cout << "depth: " << max_excess << std::endl;

        return max_excess;
    }

};