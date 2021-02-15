FLAGS=-std=c++11 -pthread -O3
BIN_FOLDER=bin
SRC_FOLDER=src
INCLUDE=$(SRC_FOLDER)/options.hpp $(SRC_FOLDER)/utils.hpp 
COMPILER=clang++

.PHONY: all 1_vector_sum 2_mmul 3_spmv 4_strings 5_bfs_dfs 6_sort

all: \
	1_vector_sum 2_mmul 3_spmv 4_strings 5_bfs_dfs 6_sort

1_vector_sum: $(INCLUDE) $(SRC_FOLDER)/1_vector_sum.cpp
	$(COMPILER) $(SRC_FOLDER)/1_vector_sum.cpp $(FLAGS) -o $(BIN_FOLDER)/1_vector_sum

2_mmul: $(INCLUDE) $(SRC_FOLDER)/2_mmul.cpp
	$(COMPILER) $(SRC_FOLDER)/2_mmul.cpp $(FLAGS) -o $(BIN_FOLDER)/2_mmul
	
3_spmv: $(INCLUDE) $(SRC_FOLDER)/3_spmv.cpp
	$(COMPILER) $(SRC_FOLDER)/3_spmv.cpp $(FLAGS) -o $(BIN_FOLDER)/3_spmv

4_strings: $(INCLUDE) $(SRC_FOLDER)/4_strings.cpp
	$(COMPILER) $(SRC_FOLDER)/4_strings.cpp $(FLAGS) -o $(BIN_FOLDER)/4_strings

5_bfs_dfs: $(INCLUDE) $(SRC_FOLDER)/5_bfs_dfs.cpp
	$(COMPILER) $(SRC_FOLDER)/5_bfs_dfs.cpp $(FLAGS) -o $(BIN_FOLDER)/5_bfs_dfs

6_sort: $(INCLUDE) $(SRC_FOLDER)/6_sort.cpp
	$(COMPILER) $(SRC_FOLDER)/6_sort.cpp $(FLAGS) -o $(BIN_FOLDER)/6_sort