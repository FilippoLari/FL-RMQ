#pragma once

#include <iostream>
#include <vector>
#include <climits>

template<typename K, typename Pos>
class SegmentTree {

    std::vector<Pos> indexes;

    size_t n;

public:

    SegmentTree() = default;

    explicit SegmentTree(const std::vector<K> &data) : n(data.size()) {
        indexes.resize(4 * n);
        rec_build(data, 0, 0, n - 1);
    }

    size_t query(const std::vector<K> &data, const size_t i, const size_t j) const {
        return rec_query(data, 0, 0, n - 1, i, j);
    }

    inline size_t size() const {
        return ((indexes.size() * sizeof(Pos)) + sizeof(n)) * CHAR_BIT;
    }

private: 

    void rec_build(const std::vector<K> &data, const size_t node,
                     const size_t left, const size_t right) {
        if (left == right) {
            indexes[node] = left;
        } else {
            int mid = (left + right) / 2;
            rec_build(data, 2 * node + 1, left, mid);
            rec_build(data, 2 * node + 2, mid + 1, right);
            indexes[node] = (data[indexes[2 * node + 1]] <= data[indexes[2 * node + 2]]) ? indexes[2 * node + 1] : indexes[2 * node + 2];
        }
    }

    size_t rec_query(const std::vector<K> &data, const size_t node,
                         const size_t left, const size_t right,
                         const size_t i, const size_t j) const {
        if (j < left || i > right) return n;
        if (i <= left && right <= j) return indexes[node];
        
        int mid = (left + right) / 2;
        int left_min_pos = rec_query(data, 2 * node + 1, left, mid, i, j);
        int right_min_pos = rec_query(data, 2 * node + 2, mid + 1, right, i, j);
        
        if (left_min_pos == n) return right_min_pos;
        if (right_min_pos == n) return left_min_pos;
        return (data[left_min_pos] <= data[right_min_pos]) ? left_min_pos : right_min_pos;
    }
};