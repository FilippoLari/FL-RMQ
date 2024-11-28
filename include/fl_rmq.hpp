// Parts of the following code are heavily adapted from: https://github.com/gvinciguerra/PGM-index
// and https://github.com/gvinciguerra/la_vector

#pragma once

#include <limits>
#include <cstdlib>
#include <climits>
#include <iostream>
#include <cassert>
#include <numeric>

#include "piecewise_linear_model.hpp"

#define SUB_EPS(x, epsilon, lo) ((x) <= (lo) + (epsilon) ? (lo) : ((x) - (epsilon)))
#define ADD_EPS(x, epsilon, hi) ((x) + (epsilon) + 2 >= (hi) ? (hi) : (x) + (epsilon) + 2)
#define SUB_DELTA(x, delta) ((x) <= delta ? 0 : ((x) - (delta)))

template<typename K, typename Range, typename Pos,
 typename Floating = float, size_t Epsilon = 64, bool Rightmost = true>
class FLRMQ {
    static_assert(std::is_integral_v<K>);
    static_assert(std::is_integral_v<Range>);
    static_assert(std::is_integral_v<Pos>);
    static_assert(std::is_floating_point_v<Floating>);
    static_assert(Epsilon > 0);

protected:

    using pla_model = OptimalPiecewiseLinearModel<Range, Pos>;

    struct Segment;

    std::vector<Segment> segments;

    std::vector<int64_t> first_segment; // todo: use the correct type
    std::vector<int64_t> deltas; // todo: use the correct type, each delta is at most nlog(n)

    size_t ranges;
    size_t n;

public:

    FLRMQ() = default;

    explicit FLRMQ(const std::vector<K> &data) : n(data.size()) {
        if(n == 0) [[unlikely]]
            return;

        std::vector<std::vector<int64_t>> st(2, std::vector<int64_t>(n));

        // initialize the PLA
        pla_model pla(Epsilon);

        segments.reserve(n/Epsilon);

        int64_t d = 0, curr = 1, prev = 0, last = n - 1, r = n;
        const auto max_p = msb(n);

        // fill the first column
        for(auto i = 0; i < n; ++i)
            st[prev][i] = i;
        
        deltas.reserve(max_p);
        deltas.push_back(d);

        for(auto j = 1; j <= max_p; ++j) {
            // rmq(i, i+(1<<j)-1) = min{rmq(i,i+(1<<(j-1))-1), rmq(i+(1<<(j-1)), i+(1<<j)-1)}
            for(auto i = 0; i + (1 << j) <= n; ++i) {
            
                if constexpr (Rightmost)
                    st[curr][i] = (data[st[prev][i]]<data[st[prev][i+(1<<(j-1))]])? st[prev][i] : st[prev][i+(1<<(j-1))];
                else
                    st[curr][i] = (data[st[prev][i]]<=data[st[prev][i+(1<<(j-1))]])? st[prev][i] : st[prev][i+(1<<(j-1))];

                // it is true O(log(n)) times
                if(i == 0 && last > st[curr][0]) [[unlikely]] {
                    d += (last - st[curr][0] + 1);
                    deltas.push_back(d);
                } else if (i == 0) [[unlikely]] {
                    deltas.push_back(d);
                }         

                last = st[curr][i];

                insert_point(r, st[curr][i] + d, pla);

                ++r;
            }
            std::swap(curr, prev);
        }

        // add the last segment
        segments.emplace_back(pla.get_segment());

        ranges = r - n;

        first_segment.reserve(max_p+1);

        // todo: can be done while building the segments
        for(auto j = 0; j <= max_p; ++j) {
            auto first = std::prev(std::upper_bound(segments.begin(), segments.end(), encode(0, ((1<<j)-1))));
            first_segment.push_back(std::distance(segments.begin(), first));
        }

        // avoid bounds checking when accessing the next diagonal at query time
        first_segment.push_back(segments.size()-1);
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

        if(len <= 2 * Epsilon + 1) { // todo: check this condition
            return find_minimum(data, i, j).second;
        }

        const auto k = msb(len);

        const auto e1 = encode(i, i + (1 << k) - 1);
        const auto e2 = encode(j - (1 << k) + 1, j);

        const auto first = std::next(segments.begin(), first_segment[k]);
        const auto last = std::next(segments.begin(), first_segment[k+1]);

        auto it_s1 = segment_for_range(e1, first, last+1);
        auto it_s2 = segment_for_range(e2, it_s1, last+1);

        auto [lo1, hi1] = reduce_interval(it_s1, e1, deltas[k], i, i + (1 << k) - 1);
        auto [lo2, hi2] = reduce_interval(it_s2, e2, deltas[k], j - (1 << k) + 1, j);

        auto [m1, p1] = find_minimum(data, lo1, hi1);
        auto [m2, p2] = find_minimum(data, lo2, hi2);

        return (m1 < m2)? p1 : p2;
    }

    /**
     * Returns the number of segments used by the PLA
     * @return the number of segments used by the PLA
     */
    inline size_t segment_count() const {
        return segments.size();
    }

    inline std::vector<Segment> get_segments() const {
        return segments;
    }

    inline size_t range_count() const {
        return ranges;
    }

    inline size_t data_count() const {
        return n;
    }

    /**
     * Returns the size in bit of the data structure
     * @return the size in bit of the data structure
     */
    inline size_t size() const {
        return ((sizeof(Segment) * segments.size()) 
                    + (deltas.size() * sizeof(int64_t))
                    + (first_segment.size() * sizeof(int64_t))
                    + sizeof(size_t) ) * CHAR_BIT;
    }

protected:

    /**
     * Returns the largest k such that 2^k <= len,
     * i.e. the most significant bit of x.
     * 
     * @param len the length of a range
     * @return the largest k such that 2^k <= len
     */
    inline uint msb(const size_t len) const { // todo: use sdsl
        const auto leading_zeroes = __builtin_clzll(len);
        return 63-leading_zeroes;
    }

    /**
     * Encodes a range whose size is a power of two in an integer
     * following the diagonal ordering.
     * 
     * @param i the left extreme of the range
     * @param j the right extreme of the range
     * @return the encoding of the range
     */
    inline size_t encode(const size_t i, const size_t j) const {
        const auto k = msb(j - i + 1); // todo: use sdsl
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
     * 2eps+1 such that the minimum of [lo,hi] is guaranteed to be inside [lo',hi'].
     * 
     * @param s the segment approximating the position of the minimum inside [lo,hi]
     * @param e the diagonal encoding of [lo,hi]
     * @param delta the cumulative correction
     * @param lo the left extreme of the interval
     * @param hi the right extreme of the interval
     * @return an interval [lo',hi'] whose size is bounded by 2eps+1, containing the minimum of [lo,hi]
     */
    inline std::pair<size_t, size_t> reduce_interval(std::vector<Segment>::const_iterator it,
                                             const size_t e, const size_t delta, const size_t lo, const size_t hi) const {
        const auto cp = SUB_DELTA((*it)(e), delta);
        const auto reduced_lo = SUB_EPS(cp, Epsilon + 1, lo);
        const auto reduced_hi = ADD_EPS(cp, Epsilon, hi);
        return std::make_pair(reduced_lo, reduced_hi);
    }

private:

    /**
     * Performs the insertion of a new point inside the PLA,
     * creating a new segment if necessary.
     * 
     * @param x the encoding of a range
     * @param y the adjusted position of a minima
     * @param pla the piecewise-linear approximation model
     */
    inline void insert_point(const Range x, const Pos y, pla_model &pla) {
        if(!pla.add_point(x, y)) {
            segments.emplace_back(pla.get_segment());
            pla.add_point(x, y);
        }
    }

    /**
     * Search for the minimum and its position inside an interval.
     * 
     * @param lo the left extreme of the interval
     * @param hi the right extreme of the interval
     * @return a pair consisting of the minimum and its position inside [lo,hi]
     */
    inline std::pair<K, size_t> find_minimum(const std::vector<K> &data, const size_t lo, const size_t hi) const {
        K min = data[lo];
        size_t idx = lo;

        for(auto i = lo; i <= hi; ++i)
            if(min > data[i]) {
                min = data[i];
                idx = i;
            }

        return std::make_pair(min, idx);
    }
};

template<typename K, typename Range, typename Pos,
 typename Floating, size_t Epsilon, bool Rightmost>
struct FLRMQ<K, Range, Pos, Floating, Epsilon, Rightmost>::Segment {
    Range range;
    Floating slope;
    Pos intercept;

    Segment() = default;

    Segment(Range range, Floating slope, Pos intercept) : range(range), slope(slope), intercept(intercept) {};

    explicit Segment(const typename OptimalPiecewiseLinearModel<Range, Pos>::CanonicalSegment &cs)
        : range(cs.get_first_x()) {
        auto [cs_slope, cs_intercept] = cs.get_floating_point_segment(range);
        if (cs_intercept > std::numeric_limits<decltype(intercept)>::max())
            throw std::overflow_error("Change the type of Segment::intercept to uint64");
        /*if (cs_intercept < 0)
            throw std::overflow_error("Unexpected intercept < 0");*/
        slope = cs_slope;
        intercept = cs_intercept;
    }

    friend inline bool operator<(const Segment &s, const Range &range) { return s.range < range; }
    friend inline bool operator<(const Range &range, const Segment &s) { return range < s.range; }
    friend inline bool operator<(const Segment &s, const Segment &t) { return s.range < t.range; }

    operator Range() const { return range; };

    /**
     * Returns the approximate position of the minimum
     * inside the interval encoded by @r.
     * 
     * @param r the encoding of a certain interval [i,j]
     * @return the approximate position of the minimum inside [i,j]
     */
    inline size_t operator()(const Range &r) const {
        size_t pos;
        if constexpr (std::is_same_v<Range, int64_t> || std::is_same_v<Range, int32_t>)
            pos = size_t(slope * double(std::make_unsigned_t<Range>(r) - range));
        else
            pos = size_t(slope * double(r - range));
        return pos + intercept;
    }

};