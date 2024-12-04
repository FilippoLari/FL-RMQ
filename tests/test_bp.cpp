#include <gtest/gtest.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <tuple>

#include "bp_utils.hpp"

using query_type = std::pair<size_t, size_t>;

std::pair<size_t, int64_t> min_excess_scan(const pasta::BitVector &bitvector,
                                            const size_t i, const size_t j) {
    int64_t min_excess = std::numeric_limits<int64_t>::max();
    size_t min_excess_idx = 0;
    int64_t excess = 0;

    for(size_t k = i; k <= j; ++k) {
        excess += (bitvector[k]) ? 1 : -1;
        min_excess = (excess <= min_excess) ? excess : min_excess;
        min_excess_idx = (excess <= min_excess) ? k : min_excess_idx;
    }

    return std::make_pair(min_excess_idx, min_excess);
}

class BitVectorTest : public ::testing::Test {
protected:

    pasta::BitVector bitvector;

    std::vector<query_type> queries;

    static constexpr int seed = 42;

    pasta::BitVector generate_bitvector(size_t size) {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<int> dis(0, 1);

        pasta::BitVector bv(size, 0);

        for (size_t i = 0; i < size; ++i) {
            bv[i] = dis(gen);
        }

        return bv;
    }

    std::vector<query_type> generate_queries(size_t num_queries, size_t size) {
    
        const std::vector<size_t> lenghts = {128, 256, 512, 1024,
                                                2048, 4096, 8192, 16384};
        std::vector<query_type> queries;
        queries.reserve(num_queries * lenghts.size());
        std::mt19937 gen(seed);

        for(const auto &length : lenghts) {
            std::uniform_int_distribution<size_t> length_dis(1, size - length + 1);

            for(auto q = 0; q < num_queries; ++q) {
                const size_t start = length_dis(gen);
                const size_t end = start + length;
                queries.emplace_back(start, end);
            }
        }

        return queries;
    }

    void SetUp() override {
        size_t size = (1 << 20);
        bitvector = generate_bitvector(size);
        queries = generate_queries(100, size);
    }
};

TEST_F(BitVectorTest, SameWordQueries) {
    for(size_t i = 0; i < 64; ++i) {
        for(size_t j = i; j < 64; ++j) {
            auto [expected_min_pos, expected_min] = min_excess_scan(bitvector, i, j);
            auto [computed_min_pos, computed_min] = min_excess(bitvector, 0, i, j);
            ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << i << ", j = " << j;
            ASSERT_EQ(expected_min, computed_min) << " Query i = " << i << ", j = " << j;
        }
    }
}

TEST_F(BitVectorTest, BetweenWordQueries) {
    for(size_t i = 0; i < 64; ++i) {
        for(size_t j = 64; j < 128; ++j) {
            auto [expected_min_pos, expected_min] = min_excess_scan(bitvector, i, j);
            auto [computed_min_pos, computed_min] = min_excess(bitvector, 0, i, j);
            ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << i << ", j = " << j;
            ASSERT_EQ(expected_min, computed_min) << " Query i = " << i << ", j = " << j;
        }
    }
}

TEST_F(BitVectorTest, GeneralQueries) {
    for(const query_type &q : queries) {
        const size_t i = q.first;
        const size_t j = q.second;
        auto [expected_min_pos, expected_min] = min_excess_scan(bitvector, i, j);
        auto [computed_min_pos, computed_min] = min_excess(bitvector, 0, i, j);
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << i << ", j = " << j;
        ASSERT_EQ(expected_min, computed_min) << " Query i = " << i << ", j = " << j;
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}