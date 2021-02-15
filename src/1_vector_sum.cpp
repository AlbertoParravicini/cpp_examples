#include <iostream>
#include <vector>
#include <chrono>
#include <numeric>
#include <thread>

#include "utils.hpp"
#include "options.hpp"

using clock_type = std::chrono::high_resolution_clock;

/////////////////////////////
/////////////////////////////

// 2020-10-18, Version 1;
template <typename T>
T sum1(std::vector<T> &a) {
    T sum = 0;
    for (auto a_i = a.begin(); a_i != a.end(); a_i++) {
        sum += *a_i;
    }
    return sum;
}

// 2020-10-18, Version 2;
template <typename T>
T sum2(std::vector<T> &a) {
    return std::accumulate(a.begin(), a.end(), (T) 0);
}

// 2020-10-18, Version 3;
#define SHIFT 128
template <typename T>
T sum3(std::vector<T> &a) {
    uint n = a.size();
    T sum[SHIFT] = { 0 };
    for (uint i = 0; i < n; i += SHIFT) {
        for (uint j = 0; j < SHIFT; j++) {
            sum[j] += (i + j) < n ? a[i + j] : (T) 0;
            sum[j] += (i + j) < n ? a[i + j] : (T) 0;
        }
    }
    T sum_final = 0;
    for (uint j = 0; j < SHIFT; j++) {
        sum_final += sum[j];
    }
    return sum_final;
}

// 2020-10-18, Version 4;
#define CORES 32

template <typename T>
void sum_inner(std::vector<T> &a, int start, int end, T &result) {
    T sum = 0;
    for (int i = start; i < end; i++) {
        sum += a[i];
    }
    result = sum;
}

template <typename T>
T sum4(std::vector<T> &a) {
    T results[CORES] = { 0 };
    std::thread threads[CORES];
    // Spawn threads;
    for (int i = 0; i < CORES; i++) {
        threads[i] = std::thread(sum_inner<T>, std::ref(a), (int) (i * (a.size() / CORES)), (int) std::min((i + 1) * (a.size() / CORES), a.size()), std::ref(results[i]));
    }
    // Join threads;
    for (int i = 0; i < CORES; i++) {
        threads[i].join();
    }
    // Final sum;
    T sum = 0;
    for (uint j = 0; j < CORES; j++) {
        sum += results[j];
    }
    return sum;
}

/////////////////////////////
/////////////////////////////

int main(int argc, char *argv[]) {
    srand(12);
    std::cout << "1 - Vector Sum;" << std::endl;

    Options options(argc, argv);
    uint n = options.N;
    uint num_iter = options.num_iter;

    std::vector<double> a(n);
    std::vector<long> times(num_iter);

    for (uint i = 0; i < num_iter; i++) {
        std::cout << "---------------- iter " << i << std::endl;
        if (i == 0 || options.reset) {
            initialize_random_array(a);
            print_array(a);
        }
        auto start = clock_type::now();
        int result = sum4(a);
        auto end = clock_type::now();

        times[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "result=" << result << std::endl;
        std::cout << "time=" << times[i] / 1000 << " ms" << std::endl;
    }

    int old_precision = std::cout.precision();
    std::cout.precision(2);
    std::cout << "----------------" << std::endl;
    std::cout << "mean execution time=" << mean<long, float>(times.data(), times.size()) / 1000 << "±" << st_dev<long, float>(times.data(), times.size()) / 1000 << " ms" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout.precision(old_precision);
}