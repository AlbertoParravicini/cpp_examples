#include <chrono>
#include <cstring>
#include <iostream>
#include <vector>
#include <limits>
#include <queue>
#include <stack>
#include <unordered_set>

#include "functions.hpp"
#include "options.hpp"
#include "utils.hpp"

using clock_type = std::chrono::high_resolution_clock;

/////////////////////////////
/////////////////////////////

// 2020-10-18, Version 1;
void bfs_1(std::vector<int> &ptr, std::vector<int> &idx, std::vector<int> &res, int start_index = 0) {
    std::queue<int> frontier;
    frontier.push(start_index);
    std::unordered_set<int> seen;
    res[start_index] = 0;
    while (frontier.size() > 0) {
        int curr_elem = frontier.front();
        frontier.pop();
        seen.insert(curr_elem);
        for (int i = ptr[curr_elem]; i < ptr[curr_elem + 1]; i++) {
            int child = idx[i];
            res[child] = std::min(res[child], res[curr_elem] + 1);
            if (seen.find(child) == seen.end()) {
                frontier.push(child);
            }
        } 
    }
}

void dfs_1(std::vector<int> &ptr, std::vector<int> &idx, std::vector<int> &res, int start_index = 0) {
    std::stack<int> stack;
    stack.push(start_index);
    std::unordered_set<int> seen;
    res[start_index] = 0;
    while (stack.size() > 0) {
        int curr_elem = stack.top();
        stack.pop();
        seen.insert(curr_elem);
        for (int i = ptr[curr_elem]; i < ptr[curr_elem + 1]; i++) {
            int child = idx[i];
            res[child] = std::min(res[child], res[curr_elem] + 1);
            if (seen.find(child) == seen.end()) {
                stack.push(child);
            }
        } 
    }
}

// dijkstra(s) {
//     S = std::set<>{s};
//     D[s] = 0;
//     P[s] = 0;
//     while S.size != N {
//         for (auto& v : S) {
//             // Neighbours of vertices in S
//             float min = INFINITY;
//             int h;
//             for (int h = ptr[*v]; h < ptr[*v + 1]; h++) {
//                 if (S.find(idx[h]) != S.end()) {
//                     min = std::min(D[v] + val[idx[h]], min);
//                 }
//             }
//             D[idx[h]] = min;
//             P[idx[h]] = *v;
//             S.push(idx[h]);
//         }
//     }
// }

// // Floyd Warshaw
// floyd_warshaw(A) {
//     float D = new[N*N];
//     int P = new[N*N];

//     // initialize D[i,j] with distances from i to j, and P[i,j] with predecessors of arcs (i,j), i.e. i; 

//     // Is going from i to j through j better than a direct jump?
//     for (int h = 0; h < N; h++){
//         for (int i = 0; i < N; i++) {
//             if (i!=h) {
//                 for (int j = 0; j < N; j++) {
//                     if (j != h) {
//                         float d_ij_new = D[i][h] + D[h][j];
//                         if (d_ij_new < D[i][j]) {
//                             D[i][j] = d_ij_new;
//                             P[i][j] = P [h][j];
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }


// ford_fulkerson(G) {
//     std::vector<int> parent(N, -1);
//     float max_flow = 0;
//     while bfs(source, sink, parent) { // Modified BFS that stores the parent edge for each vertex; edges with value < 0 are skipped; true if exist a path from source to sink
//         float path_flow = INFINITY;
//         int x = sink;
//         while (x != source) {
//             path_flow = std::min(path_flow, G[parent[x]][x]) // Min between current flow and weight of arc from parent of s to s;
//             x = parent[x];
//         }
//         max_flow += path_flow;
//         int y = sink;
//         while (y != source) {
//             int x = parent[y];
//             G[x][y] -= path_flow;
//             G[y][x] += path_flow; // If arc (x,y) now has a -1 capacity, it means that arc (y,x) now has +1 flow;
//             y = x;
//         }
//     }
//     return max_flow;
// }

/////////////////////////////
/////////////////////////////

int main(int argc, char *argv[]) {
    srand(time(0));
    std::cout << "5 - BFS & DFS;" << std::endl;

    Options options(argc, argv);
    uint n = options.N;
    uint num_iter = options.num_iter;

    std::vector<int> x;
    std::vector<int> y;
    std::vector<int> val;
    std::vector<long> times1(num_iter);
    std::vector<long> times2(num_iter);

    random_coo(x, y, val, (int)n);
    std::vector<int> ptr(n);
    std::vector<int> idx(x.size());
    coo2csr(ptr.data(), idx.data(), val.data(), x.data(), y.data(), val.data(), (int)n, (int)n, (int)x.size());

    std::vector<int> res1(n);
    std::vector<int> res2(n);

    print_graph(ptr.data(), idx.data(), n);
    print_array(val);

    int start_val = rand() % n;

    for (uint i = 0; i < num_iter; i++) {
        std::cout << "---------------- iter " << i << std::endl;
        if (i == 0 || options.reset) {
            std::fill(res1.begin(), res1.end(), std::numeric_limits<int>::max());
            std::fill(res2.begin(), res2.end(), std::numeric_limits<int>::max());
            start_val = rand() % n;
        }
        auto start = clock_type::now();
        bfs_1(ptr, idx, res1, start_val);
        auto end = clock_type::now();
        times1[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "result=" << std::endl;
        print_array(res1);
        std::cout << "time=" << times1[i] / 1000 << " ms" << std::endl;

        start = clock_type::now();
        dfs_1(ptr, idx, res2, start_val);
        end = clock_type::now();
        times2[i] = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        std::cout << "result=" << std::endl;
        print_array(res2);
        std::cout << "time=" << times2[i] / 1000 << " ms" << std::endl;

        std::cout << "errors=" << check_array_equality(res1, res2) << std::endl;
    }

    int old_precision = std::cout.precision();
    std::cout.precision(2);
    std::cout << "----------------" << std::endl;
    std::cout << "BFS mean execution time=" << mean<long, float>(times1.data(), times1.size()) / 1000 << "±" << st_dev<long, float>(times1.data(), times1.size()) / 1000 << " ms" << std::endl;
    std::cout << "DFS mean execution time=" << mean<long, float>(times2.data(), times2.size()) / 1000 << "±" << st_dev<long, float>(times2.data(), times2.size()) / 1000 << " ms" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout.precision(old_precision);
}