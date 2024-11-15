#include <iostream>
#include <numeric>
#include <fstream>
#include <sstream>

#include "libsais.h"
#include "libsais64.h"

#define LCP_EXT ".lcp"
#define BIN_EXT ".bin"
#define CSV_EXT ".csv"

/**
 * credits to: https://github.com/aboffa/Learned-Compressed-Rank-Select-TALG22/
 */
template<int First, int Last, typename Lambda>
inline void static_for_pw(Lambda const &f) {
    if constexpr (First <= Last) {
        f(std::integral_constant<size_t, First>{});
        static_for_pw<First << 1, Last>(f);
    }
}

/**
 * 
 */
template<typename T>
inline void serialize(const T *data, const size_t size, const std::string &filename) {
    std::ofstream file(filename, std::ios::binary);

    if (!file) throw std::runtime_error("Unable to serialize the data");

    file.write(reinterpret_cast<const char*>(&size), sizeof(size));

    file.write(reinterpret_cast<const char*>(data), size * sizeof(T));

    file.close();
}

/**
 * 
 */
template<typename T>
inline T* deserialize(std::ifstream &file) {
    size_t n;

    file.read(reinterpret_cast<char*>(&n), sizeof(n));

    T* data = new T[n];

    file.read(reinterpret_cast<char*>(data), n * sizeof(T));

    file.close();

    return data;
}

std::vector<int64_t> build_lcp(const std::string &filename) {

    // check if a serialized lcp already exists    
    const std::string lcp_filename = filename + LCP_EXT; 

    std::ifstream lcp_file(lcp_filename, std::ios::in | std::ios::binary);
    std::vector<int64_t> lcp_vec;
    int64_t *lcp;
    size_t n;

    if(!lcp_file) {
        std::ifstream file(filename, std::ios::in | std::ios::binary);
    
        if(!file) throw std::runtime_error("Unable to open the file");
    
        std::string s((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());

        n = s.size();

        // build the suffix array and the permuted lcp array
        uint8_t *p = reinterpret_cast<uint8_t*>(s.data());
        int64_t *sa = new int64_t[n];
        int64_t *plcp = new int64_t[n];
        lcp = new int64_t[n];

        if(libsais64(p, sa, n, 0, NULL) != 0) throw std::runtime_error("SA construction failed");

        if(libsais64_plcp(p, sa, plcp, n) != 0) throw std::runtime_error("PLCP array construction failed");
        
        if(libsais64_lcp(plcp, sa, lcp, n) != 0) throw std::runtime_error("LCP array construction failed");

        /*delete[] sa;
        delete[] p;
        delete[] plcp;*/

        // serialize the lcp array for further experiments
        serialize<int64_t>(lcp, n, lcp_filename);
    }
    else {
        //lcp = deserialize<int64_t>(lcp_file);

        lcp_file.read(reinterpret_cast<char*>(&n), sizeof(n));

        lcp = new int64_t[n];

        lcp_file.read(reinterpret_cast<char*>(lcp), n * sizeof(int64_t));

        lcp_file.close();
    }
    
    lcp_vec.assign(lcp, lcp + n);

    return lcp_vec;
}
