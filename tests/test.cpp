#include <gtest/gtest.h>
#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <tuple>

#include "fl_rmq.hpp"
#include "fl_rmq_succ.hpp"

using query_type = std::pair<size_t, size_t>;

template<typename K>
inline std::pair<K, size_t> find_minimum(const std::vector<K> &data, const size_t lo, const size_t hi) {
    K min = data[lo];
    size_t idx = lo;

    for(auto i = lo; i <= hi; ++i)
        if(data[i] < min) [[unlikely]] { // change it when testing syst version
            min = data[i];
            idx = i;
        }

    return std::make_pair(min, idx);
}

template<typename K>
std::vector<K> generate_uniform(const size_t n, const K u, const int seed = 42) {
    
    static_assert(std::is_integral<K>::value);

    std::mt19937 gen(seed); 
    
    std::uniform_int_distribution<K> dis(1, u); 

    std::vector<K> results;
    results.reserve(n);

    for(auto i = 0; i < n; ++i) {
        results.push_back(dis(gen));
    }
    
    return results;
}

std::vector<query_type> generate_queries(const size_t n, const size_t q,
                                                        const size_t r, int seed = 42) {
    if(r > n) {
        throw std::invalid_argument("Query range must be smaller than the size of the array");
    }

    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> dis(0, n - r + 1);

    std::vector<query_type> queries;
    queries.reserve(q);

    for(auto i = 0; i < q; ++i) {
        const size_t l = dis(gen);
        queries.push_back(std::make_pair(l, l+r));
    }

    return queries;
}

template<size_t Epsilon>
using FLRMQType = SuccinctFLRMQ<int, int64_t, int64_t, float, Epsilon>;

template<size_t... Epsilons>
using FLRMQTypes = std::tuple<FLRMQType<Epsilons>...>;

template<typename Tuple>
struct TupleToGTestTypes;

template<typename... Ts>
struct TupleToGTestTypes<std::tuple<Ts...>> {
    using type = ::testing::Types<Ts...>;
};

//using MyTypes = typename TupleToGTestTypes<FLRMQTypes<16, 32, 64, 128, 256, 512, 1024, 2048>>::type;

using MyTypes = typename TupleToGTestTypes<FLRMQTypes<64>>::type;

template<typename T>
class FLRMQTest : public testing::Test {
    protected:

    static std::vector<int> data;
    static std::vector<query_type> queries;

    T rmq_ds;

    static void SetUpTestSuite() {
        if(data.empty() && queries.empty()) {
            data = generate_uniform<int>(1000000, 1000000);
            std::vector<query_type> queries_50 = generate_queries(1000000, 100, 500);
            std::vector<query_type> queries_100 = generate_queries(1000000, 100, 1000);
            std::vector<query_type> queries_1000 = generate_queries(1000000, 100, 5000);
            queries.insert(queries.end(), queries_50.begin(), queries_50.end());
            queries.insert(queries.end(), queries_100.begin(), queries_100.end());
            queries.insert(queries.end(), queries_1000.begin(), queries_1000.end());
        }
    }

    FLRMQTest() : rmq_ds(data) {}
};

template<typename T>
std::vector<int> FLRMQTest<T>::data;

template<typename T>
std::vector<query_type> FLRMQTest<T>::queries;

TYPED_TEST_SUITE(FLRMQTest, MyTypes);

TYPED_TEST(FLRMQTest, QueriesWork) {
    for(const auto &q : FLRMQTest<TypeParam>::queries) {
        const auto [min, min_pos] = find_minimum(FLRMQTest<TypeParam>::data, q.first, q.second);
        const auto computed_pos = this->rmq_ds.query(q.first, q.second);
        //const auto computed_pos = this->rmq_ds.query(FLRMQTest<TypeParam>::data, q.first, q.second);
        //ASSERT_EQ(min, FLRMQTest<TypeParam>::data[computed_pos]);
        ASSERT_EQ(min_pos, computed_pos);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}