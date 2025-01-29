#include <gtest/gtest.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <tuple>

#include "sparse_table.hpp"
#include "block_decomposition.hpp"
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

typedef ::testing::Types<SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 256, 16, 256>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 512, 32, 512>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 1024, 256, 1024>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 2048, 512, 2048>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 4096, 2048, 4096>> EncodingSampledTypes;

typedef ::testing::Types<SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 256, 16>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 512, 32>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 1024, 256>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 2048, 512>,
                            SuccinctFLRMQ<int32_t, int64_t, int64_t, float, 4096, 2048>> EncodingTypes;

typedef ::testing::Types<FLRMQ<int32_t, int64_t, int64_t, float, 16, 64>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 32, 128>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 64, 256>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 128, 512>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 256, 1024>> IndexingSampledTypes;

typedef ::testing::Types<FLRMQ<int32_t, int64_t, int64_t, float, 0, 64>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 0, 128>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 0, 256>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 0, 512>,
                            FLRMQ<int32_t, int64_t, int64_t, float, 0, 1024>> IndexingTypes;

typedef ::testing::Types<BlockDecomposition<int32_t, int32_t, 10>,
                            BlockDecomposition<int32_t, int32_t, 100>,
                            BlockDecomposition<int32_t, int32_t, 1000>,
                            BlockDecomposition<int32_t, int32_t, 10000>> BlockDecompositionTypes;

typedef ::testing::Types<SparseTable<std::vector<int32_t>>> SparseTableTypes;

template <typename T>
class RMQTest : public ::testing::Test {
protected:

    static std::vector<int32_t> data;

    static std::vector<query_type> queries;

    static constexpr int seed = 42;

    static void init_data(size_t size) {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<int32_t> dis(1, size);

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
std::vector<int32_t> RMQTest<T>::data;

template<typename T>
std::vector<query_type> RMQTest<T>::queries;

template <typename T>
class EncodingSampledFLRMQTest : public RMQTest<T> {
protected:
    static void SetUpTestSuite() {
        RMQTest<T>::SetUpTestSuite();
    }
};

template <typename T>
class EncodingFLRMQTest : public RMQTest<T> {
protected:
    static void SetUpTestSuite() {
        RMQTest<T>::SetUpTestSuite();
    }
};

template <typename T>
class IndexingSampledFLRMQTest : public RMQTest<T> {
protected:
    static void SetUpTestSuite() {
        RMQTest<T>::SetUpTestSuite();
    }
};

template <typename T>
class IndexingFLRMQTest : public RMQTest<T> {
protected:
    static void SetUpTestSuite() {
        RMQTest<T>::SetUpTestSuite();
    }
};

template <typename T>
class BlockDecompositionTest : public RMQTest<T> {
protected:
    static void SetUpTestSuite() {
        RMQTest<T>::SetUpTestSuite();
    }
};

template <typename T>
class SparseTableTest : public RMQTest<T> {
protected:
    static void SetUpTestSuite() {
        RMQTest<T>::SetUpTestSuite();
    }
};

TYPED_TEST_SUITE(EncodingSampledFLRMQTest, EncodingSampledTypes);

TYPED_TEST(EncodingSampledFLRMQTest, EncodingSampledQueries) {
    TypeParam rmq_ds(EncodingFLRMQTest<TypeParam>::data);
    for(const auto &q : EncodingFLRMQTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = find_minimum(EncodingFLRMQTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(q.first, q.second);
        ASSERT_EQ(expected_min, EncodingFLRMQTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
    }
}

TYPED_TEST_SUITE(EncodingFLRMQTest, EncodingTypes);

TYPED_TEST(EncodingFLRMQTest, EncodingQueries) {
    TypeParam rmq_ds(EncodingFLRMQTest<TypeParam>::data);
    for(const auto &q : EncodingFLRMQTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = find_minimum(EncodingFLRMQTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(q.first, q.second);
        ASSERT_EQ(expected_min, EncodingFLRMQTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
    }
}

TYPED_TEST_SUITE(IndexingSampledFLRMQTest, IndexingSampledTypes);

TYPED_TEST(IndexingSampledFLRMQTest, IndexingQueries) {
    TypeParam rmq_ds(IndexingSampledFLRMQTest<TypeParam>::data);
    for(const auto &q : IndexingSampledFLRMQTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = find_minimum<int32_t, false>(IndexingSampledFLRMQTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(IndexingSampledFLRMQTest<TypeParam>::data, q.first, q.second);
        ASSERT_EQ(expected_min, IndexingSampledFLRMQTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
    }
}

TYPED_TEST_SUITE(IndexingFLRMQTest, IndexingTypes);

TYPED_TEST(IndexingFLRMQTest, IndexingQueries) {
    TypeParam rmq_ds(IndexingFLRMQTest<TypeParam>::data);
    for(const auto &q : IndexingFLRMQTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = find_minimum<int32_t, false>(IndexingFLRMQTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(IndexingFLRMQTest<TypeParam>::data, q.first, q.second);
        ASSERT_EQ(expected_min, IndexingFLRMQTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
    }
}

TYPED_TEST_SUITE(BlockDecompositionTest, BlockDecompositionTypes);

TYPED_TEST(BlockDecompositionTest, IndexingQueries) {
    TypeParam rmq_ds(BlockDecompositionTest<TypeParam>::data);
    for(const auto &q : BlockDecompositionTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = find_minimum<int32_t, true>(BlockDecompositionTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(BlockDecompositionTest<TypeParam>::data, q.first, q.second);
        ASSERT_EQ(expected_min, BlockDecompositionTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
    }
}

TYPED_TEST_SUITE(SparseTableTest, SparseTableTypes);

TYPED_TEST(SparseTableTest, IndexingQueries) {
    TypeParam rmq_ds(&SparseTableTest<TypeParam>::data);
    for(const auto &q : SparseTableTest<TypeParam>::queries) {
        const auto [expected_min, expected_min_pos] = find_minimum<int32_t, true>(SparseTableTest<TypeParam>::data, q.first, q.second);
        const auto computed_min_pos = rmq_ds.query(q.first, q.second);
        ASSERT_EQ(expected_min, SparseTableTest<TypeParam>::data[computed_min_pos]) << " Query i = " << q.first << ", j = " << q.second;
        ASSERT_EQ(expected_min_pos, computed_min_pos) << " Query i = " << q.first << ", j = " << q.second;
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}