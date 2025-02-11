// Parts of the following code are heavily adapted from: https://github.com/gvinciguerra/PGM-index
// and https://github.com/gvinciguerra/la_vector

#pragma once

#include <iostream>
#include <cstdlib>
#include <climits>
#include <cassert>
#include <numeric>
#include <limits>
#include <queue>

#include "piecewise_linear_model.hpp"

#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>

#define SUB_EPS(x, epsilon, lo) ((x) <= (lo) + (epsilon) ? (lo) : ((x) - (epsilon)))
#define ADD_EPS(x, epsilon, hi) ((x) + (epsilon) + 2 >= (hi) ? (hi) : (x) + (epsilon) + 2)
#define SUB_DELTA(x, delta) ((x) <= delta ? 0 : ((x) - (delta)))

template<typename K, typename Range, typename Pos,
 typename Floating = float, size_t Samples = 0, size_t Epsilon = 64, bool Rightmost = true>
class FLRMQ {
    static_assert(std::is_integral_v<K>);
    static_assert(std::is_integral_v<Range>);
    static_assert(std::is_integral_v<Pos>);
    static_assert(std::is_floating_point_v<Floating>);
    static_assert(Epsilon > 0);
    static_assert(Samples >= 0);

protected:

    using pla_model = OptimalPiecewiseLinearModel<Range, Pos>;

    // we avoid storing the segments as structs to decrease
    // the number of potential cache miss when searching
    // for the segment covering a given encoded range

    std::vector<Pos> intercepts;
    std::vector<Floating> slopes;

    // ranges are kept uncompressed to avoid slowing down the binary search.
    std::vector<Range> ranges;

    // they only account for O(log(n)) values hence are kept uncompressed
    std::vector<int64_t> first_segment;
    std::vector<int64_t> deltas;

    sdsl::int_vector<> samples;

    size_t n;

public:

    FLRMQ() = default;

    explicit FLRMQ(const std::vector<K> &data) : n(data.size()) {
        if(n == 0) [[unlikely]]
            return;

        std::vector<std::vector<int64_t>> st(2, std::vector<int64_t>(n));

        // initialize the PLA
        pla_model pla(Epsilon);

        slopes.reserve(n / Epsilon);
        ranges.reserve(n / Epsilon);
        intercepts.reserve(n / Epsilon);

        const auto lg_eps = sdsl::bits::hi(2 * (Epsilon + 1)); 
        const auto lg_n = sdsl::bits::hi(n);

        int64_t d = 0, curr = 1, prev = 0;

        fill_column(data, st[prev], lg_eps);

        // number of ranges originating from diagonals [0, lg_eps - 1]
        int64_t r = lg_eps * (n + 1) - (1 << lg_eps) + 1; 

        // last element of the column lg_eps - 1
        // notice that lg_eps is always >= 1
        int64_t last = st[prev][n - (1 << (lg_eps - 1)) + 1];
        
        deltas = std::vector<int64_t>(lg_n + 1, 0);

        for(auto j = lg_eps; j <= lg_n; ++j) {
            st[curr][0] = (cmp_elements(data[st[prev][0]], data[st[prev][1<<(j-1)]]))?
                                     st[prev][0] : st[prev][1<<(j-1)];

            update_deltas(st[prev], st[curr][0], lg_eps, j, d);

            for(auto i = 1; (i + (1 << j)) <= n; ++i) {
                st[curr][i] = (cmp_elements(data[st[prev][i]], data[st[prev][i+(1<<(j-1))]]))? 
                                    st[prev][i] : st[prev][i+(1<<(j-1))];

                insert_point(r, st[curr][i] + d, pla);

                ++r;
            }
            std::swap(prev, curr);
        }

        // add the last segment
        const auto [intercept, slope, range] = get_segment_params(pla);

        slopes.push_back(slope);
        ranges.push_back(range);
        intercepts.push_back(intercept);

        first_segment = std::vector<int64_t>(lg_n + 1, 0);

        // todo: can be done while building the pla
        for(auto j = lg_eps; j <= lg_n; ++j) {
            auto first =  segment_for_range(encode(0, (1 << j)-1), ranges.begin(), ranges.end());
            first_segment[j] = std::distance(ranges.begin(), first);
        }

        // avoid bounds checking when accessing the next diagonal at query time
        first_segment[lg_n] = ranges.size() - 1;

        if constexpr (Samples > 0) {
            sample(data);
        }
    }

    /**
     * Returns the position of the leftmost minimum inside the interval [i,j]
     * 
     * @param i the left extreme of the interval
     * @param j the right extreme of the interval
     * @return the position of the minimum inside [i,j] 
     */
    size_t query(const std::vector<K> &data, const size_t i, const size_t j) const {

        const auto len = j - i + 1;

        if(len <= 2 * (Epsilon + 1)) {
            return find_minimum(data, i, j).second;
        }

        const auto k = sdsl::bits::hi(len);

        const auto e1 = encode(i, i + (1 << k) - 1);
        const auto e2 = encode(j - (1 << k) + 1, j);

        const auto first = std::next(ranges.begin(), first_segment[k]);
        const auto last = std::next(ranges.begin(), first_segment[k+1]);

        const auto it_s1 = segment_for_range(e1, first, last+1);
        const auto it_s2 = segment_for_range(e2, it_s1, last+1);

        auto [lo1, hi1] = reduce_interval(std::distance(ranges.begin(), it_s1), e1,
                                             deltas[k], i, i + (1 << k) - 1);

        auto [lo2, hi2] = reduce_interval(std::distance(ranges.begin(), it_s2), e2,
                                             deltas[k], j - (1 << k) + 1, j);

        auto [m1, p1] = find_minimum(data, lo1, hi1);
        auto [m2, p2] = find_minimum(data, lo2, hi2);

        return (m1 < m2)? p1 : p2;
    }

    /**
     * Returns the number of segments used by the PLA
     * @return the number of segments used by the PLA
     */
    inline size_t segment_count() const {
        return ranges.size();
    }

    inline size_t data_count() const {
        return n;
    }

    /**
     * Returns the size in bit of the data structure
     * @return the size in bit of the data structure
     */
    inline size_t size() const {
        const auto samples_size = ((Samples > 0)? samples.bit_size() : 0);
        return ((sizeof(Range) * ranges.size() + sizeof(Pos) * intercepts.size() + sizeof(Floating) * slopes.size()) 
                    + (deltas.size() * sizeof(int64_t))
                    + (first_segment.size() * sizeof(int64_t))
                    + sizeof(size_t)) * CHAR_BIT + samples_size;
    }

protected:

    /**
     * Encodes a range whose size is a power of two in an integer
     * following the diagonal ordering.
     * 
     * @param i the left extreme of the range
     * @param j the right extreme of the range
     * @return the encoding of the range
     */
    inline size_t encode(const size_t i, const size_t j) const {
        const auto k = sdsl::bits::hi(j - i + 1);
        const auto pk = j - i + 1;
        return k * (n + 1) - pk + (i + 1);
    }

    /**
     * Search for the segment covering a given encoded range.
     * 
     * @param range the encoding of a range
     * @param first an iterator to the first segment
     * @param last an iteratore to the last segment
     * @return an iterator to the segment covering @range
     */
    template<typename It>
    inline It segment_for_range(const Range range, It first, It last) const {
        return std::prev(std::upper_bound(first, last, range));
    }

    /**
     * Given an eps-approximation for the position of the minimum inside an interval
     * [lo,hi] and its diagonal encoding, returns a smaller range [lo',hi'] of length at most
     * 2eps+2 such that the minimum of [lo,hi] is guaranteed to be inside [lo',hi'].
     * 
     * @param i the index of the segment approximating the position of the minimum inside [lo,hi]
     * @param e the diagonal encoding of [lo,hi]
     * @param delta the cumulative correction
     * @param lo the left extreme of the interval
     * @param hi the right extreme of the interval
     * @return an interval [lo',hi'] whose size is bounded by 2eps+2, containing the minimum of [lo,hi]
     */
    inline std::pair<size_t, size_t> reduce_interval(const size_t i,
                                             const size_t e, const size_t delta, const size_t lo, const size_t hi) const {
        const auto cp = SUB_DELTA(predict(i, e), delta);
        const auto reduced_lo = SUB_EPS(cp, Epsilon + 1, lo);
        const auto reduced_hi = ADD_EPS(cp, Epsilon, hi);
        return std::make_pair(reduced_lo, reduced_hi);
    }

private:

    /**
     * Perform the comparison beetween two elements a and b, based
     * on the template parameter Rightmost. Notice that a is assumed 
     * to be at the left of b.
     * 
     * @param a the first element
     * @param b the second element
     */
    inline constexpr bool compare(const K a, const K b) const {
        if constexpr (Rightmost)
            return a <= b;
        else 
            return a < b;
    }

    inline constexpr bool cmp_elements(const K a, const K b) const {
        if constexpr (Rightmost) 
            return a < b;
        else
            return a <= b;
    }

    inline void update_deltas(const std::vector<int64_t> &prev_column, const int64_t first,
                         const size_t lg_eps, const size_t j, int64_t &d) {

        const auto last = prev_column[n - (1 << (j - 1) + 1)];

        if(j != lg_eps && last > first) [[likely]] {
            d += (last - first + 1);
            deltas[j] = d;
        } else [[unlikely]] {
            deltas[j] = d;
        }
    }

    /**
     * Initialize a generic column j of a sparse table without computing all the
     * previous ones. It uses a sliding windows of size 2^(j-1) and runs in
     * linear time and space.
     * 
     * @param data the input array
     * @param column the sparse table column
     * @param j the index of the column
     */
    void fill_column(const std::vector<K> &data, std::vector<int64_t> &column,
                            const size_t j) {
        const auto w = (1 << (j - 1));

        std::deque<size_t> q;

        for(auto i = 0; i < n; ++i) {

            while(!q.empty() && q.front() + w <= i) q.pop_front();
        
            while(!q.empty() && compare(data[i], data[q.back()])) q.pop_back();

            q.push_back(i);

            if(i + 1 >= w) column[i + 1 - w] = q.front();
        }
    }

    void sample(const std::vector<K> &data) {
        const size_t samples_size = (n / Samples) + 2;

        samples = sdsl::int_vector<>(samples_size, 0);

        for(size_t i = 0; i < n; i += Samples) {

            K min = data[i];
            size_t idx = 0;

            for(size_t j = 1; i + j < n && j < Samples; ++j) {
                if(compare(data[i + j], min)) {
                    min = data[i + j];
                    idx = j;
                }
            }

            samples[i / Samples] = idx;
        }

        sdsl::util::bit_compress(samples);
    }

    /**
     * Extracts the intercept, slope and first covered range from the current
     * segment of the PLA.
     * 
     * @param pla the piecewise-linear approximation model
     */
    inline std::tuple<Pos, Floating, Range> get_segment_params(pla_model &pla) const {
        const auto segment = pla.get_segment();
        const auto range = segment.get_first_x();
        const auto [slope, intercept] = segment.get_floating_point_segment(range);

        return std::make_tuple(intercept, slope, range);
    }

    /**
     * Performs the insertion of a new point inside the PLA,
     * creating a new segment if necessary.
     * 
     * @param x the encoding of a range
     * @param y the adjusted position of a minima
     * @param pla the piecewise-linear approximation model
     * @param tmp_intercepts a the temporary vector containing all the intercepts
     *          of the segments computed so far
     */
    inline void insert_point(const Range x, const Pos y, pla_model &pla) {
        if(!pla.add_point(x, y)) {
            const auto [intercept, slope, range] = get_segment_params(pla);

            intercepts.push_back(intercept);
            slopes.push_back(slope);
            ranges.push_back(range);

            pla.add_point(x, y);
        }
    }

    inline size_t predict(const size_t i, const Range r) const {
        size_t pos;
        if constexpr (std::is_same_v<Range, int64_t> || std::is_same_v<Range, int32_t>)
            pos = size_t(slopes[i] * double(std::make_unsigned_t<Range>(r) - ranges[i]));
        else
            pos = size_t(slopes[i] * double(r - ranges[i]));
        return pos + intercepts[i];
    }

    /**
     * Search for the minimum and its position inside an interval.
     * 
     * @param lo the left extreme of the interval
     * @param hi the right extreme of the interval
     * @return a pair consisting of the minimum and its position inside [lo,hi]
     */
    inline std::pair<K, size_t> find_minimum(const std::vector<K> &data, const size_t lo, const size_t hi) const {

        if constexpr (Samples > 0) {
            
            if(hi - lo + 1 < Samples) [[unlikely]]
                return linear_scan_min(data, lo, hi);

            const size_t block_start = lo / Samples;
            const size_t block_end = hi / Samples;

            size_t block_idx = samples[block_start] + block_start * Samples;
            K min = data[block_idx];
            size_t idx = block_idx;

            for(size_t k = block_start + 1; k <= block_end; ++k) {
                block_idx = samples[k] + k * Samples;
                if(compare(data[block_idx], min)) {
                    min = data[block_idx];
                    idx = block_idx;
                }
            }
            
            if(idx >= lo && idx <= hi) [[likely]] {
                return std::make_pair(min, idx);
            } else [[unlikely]] {
                std::tie(min, idx) = linear_scan_min(data, lo, (block_start + 1) * Samples - 1);

                // todo: can be avoided
                for(size_t k = block_start + 1; k < block_end; ++k) {
                    block_idx = samples[k] + k * Samples;
                    if(compare(data[block_idx], min)) {
                        min = data[block_idx];
                        idx = block_idx;
                    }
                }
                
                const auto [last_min, last_idx] = linear_scan_min(data, block_end * Samples, hi);
                min = (compare(last_min, min)) ? last_min : min;
                idx = (compare(last_min, min)) ? last_idx : idx;

                return std::make_pair(min, idx);
            }
        } else {
            return linear_scan_min(data, lo, hi);
        }
    }

    inline std::pair<K, size_t> linear_scan_min(const std::vector<K> &data, const size_t lo, const size_t hi) const {
        K min = data[lo];
        size_t idx = lo;
        
        for(auto i = lo; i <= hi; ++i)
            if(compare(data[i], min)) {
                min = data[i];
                idx = i;
            }
        
        return std::make_pair(min, idx);
    }
};