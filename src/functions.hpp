#include <thread>

#define CORES 16

template <typename T>
inline void sum_array_inner(std::vector<T> &a, int start, int end, T &result) {
    T sum = 0;
    for (int i = start; i < end; i++) {
        sum += a[i];
    }
    result = sum;
}

template <typename T>
inline T sum_array(std::vector<T> &a) {
    T results[CORES] = {0};
    std::thread threads[CORES];
    // Spawn threads;
    for (int i = 0; i < CORES; i++) {
        threads[i] = std::thread(sum_array_inner<T>, std::ref(a), i * (a.size() / CORES), std::min((i + 1) * (a.size() / CORES), a.size()), std::ref(results[i]));
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

template <typename T>
inline void dot_product_inner(T *a, T *b, int start, int end, T *result) {
    T sum = 0;
    for (int i = start; i < end; i++) {
        sum += a[i] * b[i];
    }
    *result = sum;
}

template <typename T>
inline T dot_product(T *a, T *b, int n) {
    T results[CORES] = {0};
    std::thread threads[CORES];
    // Spawn threads;
    for (int i = 0; i < CORES; i++) {
        threads[i] = std::thread(dot_product_inner<T>, a, b, i * (n / CORES), std::min((i + 1) * (n / CORES), n), &results[i]);
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
