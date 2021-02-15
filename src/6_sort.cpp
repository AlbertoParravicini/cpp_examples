#include <algorithm>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include "options.hpp"
#include "utils.hpp"

using clock_type = std::chrono::high_resolution_clock;

/////////////////////////////
/////////////////////////////

template <typename T>
void sort_g(std::vector<T> &a) {
    std::sort(a.begin(), a.end());
}

// 2020-10-18, Version 1;
template <typename T>
bool comparator(T &lhs, T &rhs) {
    return lhs < rhs;
}

template <typename T>
void sort_1(std::vector<T> &a) {
    std::sort(a.begin(), a.end(), comparator<T>);
}

// Quicksort

// 2020-10-18, Version 1;
template <typename T>
int partition_2(std::vector<T> &a, int low, int high) {
    T pivot = a[high];
    int i = low;
    // Ensure that all values below pivot are before pivot,
    // and all values higher than pivot are after pivot;
    for (int j = low; j < high; j++) {
        if (a[j] < pivot) { // All values before j are guaranteed to be smaller than pivot
            std::swap(a[i], a[j]);
            i++;
        }
    }
    // Put the pivot in the right place;
    std::swap(a[high], a[i]);
    return i;
}

template <typename T>
void sort_2_inner(std::vector<T> &a, int low, int high) {
    if (low < high) {
        // Place the pivot in the right sorted place;
        int p = partition_2(a, low, high);
        // Now apply to the remaining part of the arrays;
        sort_2_inner(a, low, p - 1);
        sort_2_inner(a, p + 1, high);
    }
}

template <typename T>
void sort_2(std::vector<T> &a) {
    sort_2_inner(a, 0, a.size() - 1);
}

// 2020-10-18, Version 2;
template <typename T>
int partition_3(std::vector<T> &a, int low, int high) {
    int m = (high + low) / 2;
    T pivot = a[m];
    // if (a[m] < a[low]) std::swap(a[low], a[m]); // If low is bigger than middle, put low in the center;
    // if (a[high] < a[low]) std::swap(a[high], a[low]); // If low is bigger than high, swap them;
    // if (a[m] < a[high]) std::swap(a[high], a[m]);// If high is bigger than middle, put high in the center;
    // pivot = a[high]; // The median is now stored at high;

    // Initialize with -+1 as we have a do-while;
    int i = low - 1;
    int j = high + 1;
    // Ensure that all values below pivot are before pivot,
    // and all values higher than pivot are after pivot;
    while (true) {
        // Do-while to skip after the current values;
        do {
            i++;
        } while (a[i] < pivot);
        do {
            j--;
        } while (a[j] > pivot);
        if (i >= j)
            return j;
        std::swap(a[j], a[i]);
    }
}

template <typename T>
void sort_3_inner(std::vector<T> &a, int low, int high) {
    if (low < high) {
        // The pivot is not necessarily placed at p here!
        int p = partition_3(a, low, high);
        // Now apply to the remaining part of the arrays.
        // Tail-recursion on the larger partition;
        if (p - low < high - p - 1) {
            sort_3_inner(a, low, p);
            sort_3_inner(a, p + 1, high);
        } else {
            sort_3_inner(a, p + 1, high);
            sort_3_inner(a, low, p);
        }
    }
}

template <typename T>
void sort_3(std::vector<T> &a) {
    sort_3_inner(a, 0, a.size() - 1);
}

// Insertion Sort;

// 2020-10-18, Version 4;
template <typename T>
void sort_4(T *a, int start, int end) {
    for (int i = start + 1; i < end; i++) {
        int j = i;
        while ((j > 0) && (a[j] < a[j - 1])) {
            std::swap(a[j], a[j - 1]);
            j--;
        } 
    }
}

template <typename T>
void sort_4(std::vector<T> &a) {
    sort_4(a.data(), 0, a.size());
}

// 2020-10-18, Version 5;
template <typename T, int threshold>
void sort_5_inner(std::vector<T> &a, int low, int high) {
    if (low < high) {
        // The pivot is not necessarily placed at p here!
        int p = partition_3(a, low, high);
        // Now apply to the remaining part of the arrays.
        if (p - low > threshold) {
            sort_5_inner<T, threshold>(a, low, p);
        } else {
            sort_4(a.data(), low, p + 1);
        }
        if (high - p - 1 > threshold) {
            sort_5_inner<T, threshold>(a, p + 1, high);
        } else {
            sort_4(a.data(), p + 1, high + 1);
        }
    }
}

template <typename T>
void sort_5(std::vector<T> &a) {
    sort_5_inner<T, 4000>(a, 0, a.size() - 1);
}

// 2020-10-18, Version 6;

#define THREADS 16
#define THRESHOLD_2 4000
template <typename T>
void sort_6(std::vector<T> &a) {
    std::thread threads[THREADS];
    for (int i = 0; i < THREADS; i++) {
        int low = std::min(i * ((a.size() + THREADS - 1) / THREADS), a.size() - 1);
        int high =  std::min((i + 1) * ((a.size() + THREADS - 1) / THREADS) - 1, a.size() - 1);
        threads[i] = std::thread(sort_5_inner<T, 0>, std::ref(a), low, high);
        // std::cout << i << " " << low << " " << high << std::endl;
        // print_array(a);
        // sort_3_inner(a, low, high);
    }
    for (int i = 0; i < THREADS; i++) {
        threads[i].join();
    }
    // Final sort;
    sort_5(a);
}

// 2020-10-18, Version 7;

// #define THREADS 16
// template <typename T>
// void sort_7(std::vector<T> &a) {
//     std::thread threads[THREADS];
//     for (int i = 0; i < THREADS; i++) {
//         int low = std::min(i * ((a.size() + THREADS - 1) / THREADS), a.size() - 1);
//         int high =  std::min((i + 1) * ((a.size() + THREADS - 1) / THREADS) - 1, a.size() - 1);
//         threads[i] = std::thread(std::sort<T>, a.data() + low, a.data() + high);
//     }
//     for (int i = 0; i < THREADS; i++) {
//         threads[i].join();
//     }
//     // Final sort;
//     std::sort(a.begin(), a.end());
// }

/////////////////////////////
/////////////////////////////

int main(int argc, char *argv[]) {
    srand(12);
    std::cout << "1 - Vector Sum;" << std::endl;

    Options options(argc, argv);
    uint n = options.N;
    uint num_iter = options.num_iter;

    std::vector<int> a(n);
    std::vector<int> a_g = a;
    std::vector<long> times1(num_iter);
    std::vector<long> times2(num_iter);

    for (uint i = 0; i < num_iter; i++) {
        std::cout << "---------------- iter " << i << std::endl;
        if (i == 0 || options.reset) {
            initialize_random_array(a);
            a_g = a;
            print_array(a);
        }
        auto start = clock_type::now();
        sort_g(a_g);
        auto end = clock_type::now();
        times1[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        start = clock_type::now();
        sort_6(a);
        end = clock_type::now();

        times2[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "result=" << check_array_equality(a, a_g) << std::endl;
        print_array(a_g);
        print_array(a);
        std::cout << "time, gold=" << times1[i] / 1000 << " ms" << std::endl;
        std::cout << "time=" << times2[i] / 1000 << " ms" << std::endl;
    }

    int old_precision = std::cout.precision();
    std::cout.precision(2);
    std::cout << "----------------" << std::endl;
    std::cout << "mean execution time, gold=" << mean<long, float>(times1.data(), times1.size()) / 1000 << "±" << st_dev<long, float>(times1.data(), times1.size()) / 1000 << " ms" << std::endl;
    std::cout << "mean execution time=" << mean<long, float>(times2.data(), times2.size()) / 1000 << "±" << st_dev<long, float>(times2.data(), times2.size()) / 1000 << " ms" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout.precision(old_precision);
}