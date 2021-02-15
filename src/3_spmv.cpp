#include <chrono>
#include <cstring>
#include <iostream>
#include <vector>

#include "functions.hpp"
#include "options.hpp"
#include "utils.hpp"

using clock_type = std::chrono::high_resolution_clock;

/////////////////////////////
/////////////////////////////

// a is N x K;
// b is K x M;
// c is N x M;

// 2020-10-18, Version 1;
template <typename I, typename T>
void spmv_csr1(std::vector<I> &ptr, std::vector<I> &idx, std::vector<T> &val, std::vector<T> &vec, std::vector<T> &res) {
    for (int i = 0; i < ptr.size(); i++) {
        T result = 0;
        for (int j = ptr[i]; j < ptr[i + 1]; j++) {
            result += vec[idx[j]] * val[j];
        }
        res[i] = result;
    }
}

// 2020-10-18, Version 2;
#define PARTITIONS 16
#define INNER_PARTITIONS 2
template <typename I, typename T>
void spmv_inner(std::vector<I> &ptr, std::vector<I> &idx, std::vector<T> &val, std::vector<T> &vec, std::vector<T> &res, int start, int end) {
    for (int i = start; i < end; i++) {
        T result = 0;
        for (int j = ptr[i]; j < ptr[i + 1]; j++) {
            result += vec[idx[j]] * val[j];
        }
        res[i] = result;
    }
}

// template <typename I, typename T>
// void spmv_inner(std::vector<I> &ptr, std::vector<I> &idx, std::vector<T> &val, std::vector<T> &vec, std::vector<T> &res, int start, int end) {
//     int N = val.size();
//     for (int i = start; i < end; i++) {
//         T result = 0;
//         for (int j = ptr[i]; j < ptr[i + 1]; j += INNER_PARTITIONS) {
//             for (int k = 0; k < INNER_PARTITIONS; k++) {
//                 result += (j + k < N) ? (vec[idx[j + k]] * val[j + k]) : (T) 0;
//             }
//         }
//         res[i] = result;
//     }
// }

template <typename I, typename T>
void spmv_csr2(std::vector<I> &ptr, std::vector<I> &idx, std::vector<T> &val, std::vector<T> &vec, std::vector<T> &res) {
    std::thread threads[PARTITIONS];
    int N = ptr.size();
    for (int i = 0; i < PARTITIONS; i++) {
        int start = i * (N + PARTITIONS - 1) / PARTITIONS;
        int end = (i + 1) * (N + PARTITIONS - 1) / PARTITIONS;
        threads[i] = std::thread(spmv_inner<I, T>, std::ref(ptr), std::ref(idx), std::ref(val), std::ref(vec), std::ref(res), start, end);
    }
    for (int i = 0; i < PARTITIONS; i++) {
        threads[i].join();
    }
}

// 2020-10-18, Version 1;
template <typename I, typename T>
void spmv_coo1(std::vector<I> &x, std::vector<I> &y, std::vector<T> &val, std::vector<T> &vec, std::vector<T> &res) {
    I x_old = 0;
    T result = 0;
    for (int i = 0; i < x.size(); i++) {
        I x_curr = x[i];
        T res_tmp = val[i] * vec[y[i]];
        if (x_curr != x_old) {
            res[x_old] = result;
            x_old = x_curr;
            result = 0;
        }
        result += res_tmp;
    }
    // Final value;
    res[x_old] = result;
}

// 2021-02-08, Version 2;
template <typename I, typename T>
void spmv_coo2(std::vector<I> &x, std::vector<I> &y, std::vector<T> &val, std::vector<T> &vec, std::vector<T> &res) {
    for (int i = 0; i < x.size(); i++) {
        res[x[i]] += val[i] * vec[y[i]];
    }
}

/////////////////////////////
/////////////////////////////

int main(int argc, char *argv[]) {
    srand(12);
    std::cout << "3 - SpMV;" << std::endl;

    Options options(argc, argv);
    int n = options.N;
    int nnz = 0;
    int num_iter = options.num_iter;

    std::vector<int> x;
    std::vector<int> y;
    std::vector<float> val;
    std::vector<long> times1(num_iter);
    std::vector<long> times2(num_iter);

    // random_coo(x, y, val, (int)n);
    readMtx<int, float>("data/coAuthorsCiteseer.mtx", &x, &y, &val, &n, &n, &nnz, 1, false, options.debug, false, true);
    
    std::vector<float> vec(n);
    std::vector<float> res(n);

    std::vector<int> ptr(n);
    std::vector<int> idx(x.size());
    coo2csr(ptr.data(), idx.data(), val.data(), x.data(), y.data(), val.data(), (int)n, (int)n, (int)x.size());

    print_graph(ptr.data(), idx.data(), n);
    print_array(val);

    for (uint i = 0; i < num_iter; i++) {
        std::cout << "---------------- iter " << i << std::endl;
        if (i == 0 || options.reset) {
            initialize_random_array(vec);
            print_array(vec);
        }
        auto start = clock_type::now();
        spmv_csr1(ptr, idx, val, vec, res);
        // spmv_coo1(x, y, val, vec, res);
        auto end = clock_type::now();
        times1[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "result=" << std::endl;
        print_array(res);
        std::cout << "time=" << (float) times1[i] / 1000 << " ms" << std::endl;

        start = clock_type::now();
        spmv_coo2(x, y, val, vec, res);
        end = clock_type::now();
        times2[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "result=" << std::endl;
        print_array(res);
        std::cout << "time=" << (float) times2[i] / 1000 << " ms" << std::endl;
    }

    int old_precision = std::cout.precision();
    std::cout.precision(2);
    std::cout << "----------------" << std::endl;
    std::cout << "CSR mean execution time=" << mean<long, float>(times1.data(), times1.size()) / 1000 << "±" << st_dev<long, float>(times1.data(), times1.size()) / 1000 << " ms" << std::endl;
    std::cout << "COO mean execution time=" << mean<long, float>(times2.data(), times2.size()) / 1000 << "±" << st_dev<long, float>(times2.data(), times2.size()) / 1000 << " ms" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout.precision(old_precision);
}