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
template <typename T>
void mmul1(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c, int N, int K, int M) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            for (int q = 0; q < K; q++) {
                c[i * M + j] += a[i * K + q] * b[q * M + j];
            }
        }
    }
}

// 2020-10-18, Version 2;
template <typename T>
void mmul2(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c, int N, int K, int M) {

    // Transpose b, now it is M x K;
    std::vector<T> b_t(K * M);

    for (int j = 0; j < M; j++) {
        for (int q = 0; q < K; q++) {
            b_t[j * K + q] = b[q * M + j];
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            T sum = 0;
            for (int q = 0; q < K; q++) {
                sum += a[i * K + q] * b_t[j * K + q];
            }
            c[i * M + j] = sum;
        }
    }
}

// 2020-10-18, Version 3;
template <typename T>
void mmul3(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c, int N, int K, int M) {

    // Transpose b, now it is M x K;
    std::vector<T> b_t(K * M);

    for (int j = 0; j < M; j++) {
        for (int q = 0; q < K; q++) {
            b_t[j * K + q] = b[q * M + j];
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            c[i * M + j] = dot_product(a.data() + i * K, b_t.data() + j * K, K);
        }
    }
}

// 2020-10-18, Version 4;
#define PARTITIONS 6

template <typename T>
void mmul4(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c, int N, int K, int M) {

    // Transpose b, now it is M x K;
    std::vector<T> b_t(K * M);

    for (int j = 0; j < M; j++) {
        for (int q = 0; q < K; q++) {
            b_t[j * K + q] = b[q * M + j];
        }
    }

    for (int i = 0; i < PARTITIONS; i++) {
        for (int j = 0; j < PARTITIONS; j++) {
            int i_start = i * ((N + PARTITIONS - 1) / PARTITIONS);
            int i_end = (i + 1) * ((N + PARTITIONS - 1) / PARTITIONS);
            int j_start = j * ((M + PARTITIONS - 1) / PARTITIONS);
            int j_end = (j + 1) * ((M + PARTITIONS - 1) / PARTITIONS);

            for (int i = i_start; i < std::min(i_end, N); i++) {
                for (int j = j_start; j < std::min(j_end, M); j++) {
                    T sum = 0;
                    for (int q = 0; q < K; q++) {
                        sum += a[i * K + q] * b_t[j * K + q];
                    }
                    c[i * M + j] = sum;
                }
            }
        }
    }
}

// 2020-10-18, Version 5;

template <typename T>
void mmul_inner(T *a, T *b, T *c, int N_start, int N_end, int M_start, int M_end, int K, int M) {
    for (int i = N_start; i < N_end; i++) {
        for (int j = M_start; j < M_end; j++) {
            T sum = 0;
            for (int q = 0; q < K; q++) {
                sum += a[i * K + q] * b[j * K + q];
            }
            c[i * M + j] = sum;
        }
    }
}

template <typename T>
void mmul5(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c, int N, int K, int M) {

    // Transpose b, now it is M x K;
    std::vector<T> b_t(K * M);

    for (int j = 0; j < M; j++) {
        for (int q = 0; q < K; q++) {
            b_t[j * K + q] = b[q * M + j];
        }
    }

    std::thread threads[PARTITIONS * PARTITIONS];
    for (int i = 0; i < PARTITIONS; i++) {
        for (int j = 0; j < PARTITIONS; j++) {
            int i_start = i * ((N + PARTITIONS - 1) / PARTITIONS);
            int i_end = (i + 1) * ((N + PARTITIONS - 1) / PARTITIONS);
            int j_start = j * ((M + PARTITIONS - 1) / PARTITIONS);
            int j_end = (j + 1) * ((M + PARTITIONS - 1) / PARTITIONS);

            i_end = std::min(i_end, N);
            j_end = std::min(j_end, M);

            threads[i * PARTITIONS + j] = std::thread(mmul_inner<T>, a.data(), b_t.data(), c.data(), i_start, i_end, j_start, j_end, K, M);
        }
    }

    for (int i = 0; i < PARTITIONS * PARTITIONS; i++) {
        threads[i].join();
    }
}

// 2020-10-20, Version 6;
#define BLOCK_SIZE 16
template <typename T>
void mmul6(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c, int N, int K, int M) {
    for (int k = 0; k < K; k += BLOCK_SIZE) {
        for (int j = 0; j < M; j += BLOCK_SIZE) {
            for (int i = 0; i < N; i++) {
                for (int jj = j; jj < std::min(j + BLOCK_SIZE, M); jj++) {
                    for (int kk = k; kk < std::min(k + BLOCK_SIZE, K); kk++) {
                        c[i * M + jj] += a[i * K + kk] * b[kk * K + jj];
                    }
                }
            }
        }
    }
}

template <typename T>
void mmul7(std::vector<T> &a, std::vector<T> &b, std::vector<T> &c, int N, int K, int M) {
    
    // Transpose b, now it is M x K;
    std::vector<T> b_t(K * M);

    for (int j = 0; j < M; j++) {
        for (int q = 0; q < K; q++) {
            b_t[j * K + q] = b[q * M + j];
        }
    }

    for (int k = 0; k < K; k += BLOCK_SIZE) {
        for (int j = 0; j < M; j += BLOCK_SIZE) {
            for (int i = 0; i < N; i++) {
                for (int jj = j; jj < std::min(j + BLOCK_SIZE, M); jj++) {
                    for (int kk = k; kk < std::min(k + BLOCK_SIZE, K); kk++) {
                        c[i * M + jj] += a[i * K + kk] * b_t[jj * M + kk];
                    }
                }
            }
        }
    }
}


/////////////////////////////
/////////////////////////////

int main(int argc, char *argv[]) {
    srand(12);
    std::cout << "2 - Matrix Multiplication;" << std::endl;

    Options options(argc, argv);
    uint n = options.N;
    uint num_iter = options.num_iter;

    std::vector<int> a(n * n);
    std::vector<int> b(n * n);
    std::vector<int> c(n * n);
    std::vector<long> times(num_iter);

    for (uint i = 0; i < num_iter; i++) {
        std::cout << "---------------- iter " << i << std::endl;
        if (i == 0 || options.reset) {
            initialize_random_array(a);
            print_matrix(a, n, n);
            initialize_random_array(b);
            print_matrix(b, n, n);
        }
        memset((void *)c.data(), 0, n * n * sizeof(int));
        auto start = clock_type::now();
        mmul7(a, b, c, n, n, n);
        auto end = clock_type::now();

        times[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "result=" << std::endl;
        print_matrix(c, n, n);
        std::cout << "time=" << times[i] / 1000 << " ms" << std::endl;
    }

    int old_precision = std::cout.precision();
    std::cout.precision(2);
    std::cout << "----------------" << std::endl;
    std::cout << "mean execution time=" << mean<long, float>(times.data(), times.size()) / 1000 << "±" << st_dev<long, float>(times.data(), times.size()) / 1000 << " ms" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout.precision(old_precision);
}