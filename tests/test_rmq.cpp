#include <gtest/gtest.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <tuple>

#include "fl_rmq_succ.hpp"
#include "fl_rmq.hpp"

using query_type = std::pair<size_t, size_t>;

template<typename K, bool Leftmost = true>
inline bool compare(const K &a, const K &b) {
    if constexpr (Leftmost) 
        return a < b;
    else
        return a <= b;
}

template<typename K, bool Leftmost = true>
inline std::pair<K, size_t> find_minimum(const std::vector<K> &data, const size_t lo, const size_t hi) {
    K min = data[lo];
    size_t idx = lo;

    for(auto i = lo; i <= hi; ++i)
        if(compare<K, Leftmost>(data[i], min)) {
            min = data[i];
            idx = i;
        }

    return std::make_pair(min, idx);
}

typedef ::testing::Types<SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 256, 16>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 512, 32>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 1024, 256>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 2048, 512>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 4096, 2048>> NonSystTypes;

typedef ::testing::Types<FLRMQ<int32_t, int64_t, int64_t, float, 64>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 128>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 256>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 512>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 1024>> SystTypes;

template <typename T>
class FLRMQTest : public ::testing::Test {
protected:

    static std::vector<int32_t> data;

    static std::vector<query_type> queries;

    static constexpr int seed = 42;

    static void init_data(size_t size) {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<int32_t> dis(1, 1e9);

        for(size_t i = 0; i < size; ++i)
            data.push_back(dis(gen));
    }

    static void init_queries(size_t num_queries, size_t size) {
        const std::vector<size_t> lenghts = {200, 600, 1200, 
                                                6000, 8000, 10000};
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
    }

    static void SetUpTestSuite() {
        if(queries.empty()) {
            size_t size = 1e6;
            init_data(size);
            init_queries(100, size);
        }
    }
};

template<typename T>
std::vector<int32_t> FLRMQTest<T>::data;

template<typename T>
std::vector<query_type> FLRMQTest<T>::queries;

template <typename T>
class NonSystFLRMQTest : public FLRMQTest<T> {
protected:
    static void SetUpTestSuite() {
        FLRMQTest<T>::SetUpTestSuite();
    }
};

template <typename T>
class SystFLRMQTest : public FLRMQTest<T> {
protected:
    static void SetUpTestSuite() {
        FLRMQTest<T>::SetUpTestSuite();
    }
};

TYPED_TEST_SUITE(NonSystFLRMQTest, NonSystTypes);

TYPED_TEST(NonSystFLRMQTest, NonSystQueries) {
    TypeParam rmq_ds(NonSystFLRMQTest<TypeParam>::data);
    for(const auto &q : NonSystFLRMQTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = find_minimum(NonSystFLRMQTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(q.first, q.second);
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min, NonSystFLRMQTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
    }
}

/*
TYPED_TEST_SUITE(SystFLRMQTest, SystTypes);

TYPED_TEST(SystFLRMQTest, SystQueries) {
    TypeParam rmq_ds(SystFLRMQTest<TypeParam>::data);
    for(const auto &q : SystFLRMQTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = find_minimum<int32_t, false>(SystFLRMQTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(SystFLRMQTest<TypeParam>::data, q.first, q.second);
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min, SystFLRMQTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
    }
}*/

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}