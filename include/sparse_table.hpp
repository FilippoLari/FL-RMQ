#pragma once

#include <limits>
#include <cstdlib>
#include <climits>
#include <numeric>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/bits.hpp>
#include <sdsl/util.hpp>

template<class Container>
class SparseTable
{
    std::vector<sdsl::int_vector<>> table;
    const Container *data;
    
public:
    
    SparseTable() = default;

    explicit SparseTable(const Container *data) : data(data) {
        if(data == nullptr) [[unlikely]]
            return;

        const size_t n = data->size();

        const auto max_p = sdsl::bits::hi(n);
        table.resize(max_p);

        for (auto i = 0; i < max_p; ++i) {
            table[i] = sdsl::int_vector<>(n - (1ULL << (i + 1)) + 1, 0, i + 1);
        }

        for (auto i = 0; i < n - 1; ++i) {
            if((*data)[i+1] < (*data)[i])
                table[0][i] = 1;
        }

        for (auto i = 1; i < max_p; ++i) {
            for (auto j = 0; j < table[i].size(); ++j) {
                table[i][j] = ((*data)[j+table[i-1][j]] <= (*data)[j+(1ULL<<i)+table[i-1][j+(1ULL<<i)]]) 
                                ? table[i-1][j] : (1ULL<<i)+table[i-1][j+(1ULL<<i)];           
            }
        }

        for(auto &row : table)
            sdsl::util::bit_compress(row);
    }

    inline size_t query(const size_t i, const size_t j) const {
        if (i == j) [[unlikely]]
            return i;
        if(i + 1 == j) [[unlikely]]
            return ((*data)[i] < (*data)[j]) ? i : j;
        const auto k = sdsl::bits::hi(j - i + 1);
        const auto jj = j - (1ULL << k) + 1;
        return ((*data)[i+table[k-1][i]] <= (*data)[jj+table[k-1][jj]]) ? i+table[k-1][i] : jj+table[k-1][jj];
    }

    inline size_t size() const {
        auto table_size = 0;
        for(const auto &row : table)
            table_size += row.bit_size();
        return table_size + (sizeof(data) * CHAR_BIT);
    }
};
