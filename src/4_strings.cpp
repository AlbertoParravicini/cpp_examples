#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <unordered_map>

#include "utils.hpp"
#include "options.hpp"

using clock_type = std::chrono::high_resolution_clock;

/////////////////////////////
/////////////////////////////

// 2020-10-18, Version 1;
std::string reverse(std::string s) {
    std::string s_r;
    for (auto c = s.end(); c != s.begin(); ) {
        c--;
        s_r += *c;
    }
    return s_r;
}

std::string reverse2(std::string s) {
    std::string s_r;
    for (auto c = s.rbegin(); c != s.rend(); c++) {
        s_r += *c;
    }
    return s_r;
}

std::string reverse3(std::string s) {
    return std::string(s.rbegin(), s.rend());
}

bool palyndrome(std::string s) {
    for (int i = 0; i < s.length() / 2; i++) {
        if (s[i] != s[s.length() - 1 - i]) return false;
    }
    return true;
}

bool palyndrome_exists(std::string s) {
    std::unordered_map<char, int> map;
    for (auto c : s) {
        if (map.find(c) == map.end()) {
            map[c] = 1;
        } else {
            map[c]++;
        }
    }
    bool odd_found = false;
    for (auto it : map) {
        if (it.second % 2) {
            if (!odd_found) odd_found = true;
            else return false;
        }
    }
    return true;
}


/////////////////////////////
/////////////////////////////

int main(int argc, char *argv[]) {
    srand(12);
    std::cout << "1 - String operations;" << std::endl;

    Options options(argc, argv);
    uint n = options.N;
    uint num_iter = options.num_iter;

    std::string s;
    std::vector<long> times(num_iter);

    for (uint i = 0; i < num_iter; i++) {
        std::cout << "---------------- iter " << i << std::endl;
        if (i == 0 || options.reset) {
            s = random_string(n);
            s += "ab" + reverse(s);
            // s = "aabbcddc";
            print_string(s);
        }
        auto start = clock_type::now();
        bool result = palyndrome_exists(s);
        auto end = clock_type::now();

        times[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "result = " << result << std::endl;
        std::cout << "time=" << times[i] / 1000 << " ms" << std::endl;
    }

    int old_precision = std::cout.precision();
    std::cout.precision(2);
    std::cout << "----------------" << std::endl;
    std::cout << "mean execution time=" << mean<long, float>(times.data(), times.size()) / 1000 << "±" << st_dev<long, float>(times.data(), times.size()) / 1000 << " ms" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout.precision(old_precision);
}