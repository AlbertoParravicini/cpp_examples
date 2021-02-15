#pragma once
#include <iostream>
#include <random>
#include <tuple>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include "mmio.hpp"

#define int_type int

template <typename T>
inline void initialize_random_array(T *array, int n, T min_val = 0, T max_val = 10, bool sum_to_one = true, bool norm_one = false) {
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_real_distribution<double> dist(min_val, max_val);
    for (int i = 0; i < n; i++) {
        array[i] = (T)dist(engine);
    }

    if (sum_to_one) {
        T sum = 0;
        for (int i = 0; i < n; i++) {
            sum += array[i];
        }
        for (int i = 0; i < n; i++) {
            array[i] /= sum;
        }
    } else if (norm_one) {
        T sum = 0;
        for (int i = 0; i < n; i++) {
            sum += array[i] * array[i];
        }
        for (int i = 0; i < n; i++) {
            array[i] /= sqrt((double)sum);
        }
    }
}

template <typename T>
inline void initialize_random_array(std::vector<T> &array, T min_val = 0, T max_val = 10, bool sum_to_one = true, bool norm_one = false) {
    initialize_random_array(array.data(), array.size(), min_val, max_val, sum_to_one, norm_one);
}

inline void initialize_random_array(int *array, int n, int min_val = 0, int max_val = 10) {
    for (int i = 0; i < n; i++) {
        array[i] = (rand() % (max_val - min_val)) + min_val;
    }
}

inline void initialize_random_array(std::vector<int> &array, int min_val = 0, int max_val = 10) {
    initialize_random_array(array.data(), array.size(), min_val, max_val);
}

template <typename T>
inline void print_array(T *array, int n, bool indexed = false, int max_print = 20) {
    if (indexed) {
        std::cout << "[" << std::endl;
        if (n <= max_print) {
            for (int i = 0; i < n; i++) {
                std::cout << i << ") " << array[i] << ", " << std::endl;
            }
        } else {
            int print_size = max_print / 2;
            for (int i = 0; i < print_size; i++) {
                std::cout << i << ") " << array[i] << ", " << std::endl;
            }
            std::cout << "(...),";
            for (int i = n - print_size - 1; i < n; i++) {
                std::cout << i << ") " << array[i] << ", " << std::endl;
            }
        }
        std::cout << "]" << std::endl;
    } else {
        std::cout << "[";
        if (n <= max_print) {
            for (int i = 0; i < n; i++) {
                std::cout << array[i];
                if (i < n - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
        } else {
            int print_size = max_print / 2;
            for (int i = 0; i < print_size; i++) {
                std::cout << array[i] << ", ";
            }
            std::cout << "..., ";
            for (int i = n - print_size - 1; i < n; i++) {
                std::cout << array[i];
                if (i < n - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
        }
    }
}

template <typename T>
inline void print_array(std::vector<T> &vector, bool indexed = false, int max_print = 20) {
    print_array(vector.data(), vector.size(), indexed, max_print);
}

template <typename T>
inline void print_matrix(std::vector<T> &vector, uint n, uint m, bool indexed = false, int max_print = 20) {
    // Row-wise print;
    std::cout << "[" << std::endl;
    int print_size = max_print / 2;
    if (n <= max_print) {
        for (uint i = 0; i < n; i++) {
            if (indexed) {
                std::cout << "---- row=" << i << std::endl;
            }
            print_array(vector.data() + m * i, m, indexed, max_print);
        }
    } else {
        for (uint i = 0; i < print_size; i++) {
            if (indexed) {
                std::cout << "---- row=" << i << std::endl;
            }
            print_array(vector.data() + m * i, m, indexed, max_print);
        }
        std::cout << "(...)," << std::endl;
        for (uint i = n - print_size - 1; i < n; i++) {
            if (indexed) {
                std::cout << "---- row=" << i << std::endl;
            }
            print_array(vector.data() + m * i, m, indexed, max_print);
        }
    }
    std::cout << "]" << std::endl;
}

template <typename I, typename O>
inline O mean(I *a, int n, int skip = 0) {
    O sum = 0;
    int fixed_size = n - skip;
    if (fixed_size <= 0)
        return (O)0;
    for (int i = skip; i < n; i++) {
        sum += (O)a[i];
    }
    return sum / fixed_size;
}

template <typename I, typename O>
inline O st_dev(I *a, int n, int skip = 0) {
    O sum = 0;
    O sum_sq = 0;
    int fixed_size = n - skip;
    if (fixed_size <= 0)
        return (O)0;
    for (int i = skip; i < n; i++) {
        sum += (O)a[i];
        sum_sq += (O)a[i] * (O)a[i];
    }
    O diff = sum_sq - (sum * sum / fixed_size);
    return diff >= 0 ? std::sqrt(diff / fixed_size) : (O)0;
}

template <typename I, typename T>
inline void random_coo(std::vector<I> &x, std::vector<I> &y, std::vector<T> &val, I N, I degree_max = 10, T min_value = 0, T max_value = 10) {

    std::random_device random;
    std::mt19937 engine(random());
    std::uniform_real_distribution<double> dist(min_value, max_value);

    for (int i = 0; i < N; i++) {
        std::set<I> edges;
        while (edges.size() < rand() % (degree_max + 1)) {
            edges.insert((I) rand() % N);
        }
        int j = 0;
        for (auto iter = edges.begin(); iter != edges.end(); iter++) {
            x.push_back(i);
            y.push_back(*iter);
            val.push_back((T) dist(engine));
        }
    }
}

template <typename I>
inline void random_coo(std::vector<I> &x, std::vector<I> &y, std::vector<int> &val, I N, I degree_max = 10, int min_value = 0, int max_value = 10) {
    for (int i = 0; i < N; i++) {
        std::set<I> edges;
        while (edges.size() < rand() % (degree_max + 1)) {
            edges.insert((I) rand() % N);
        }
        int j = 0;
        for (auto iter = edges.begin(); iter != edges.end(); iter++) {
            x.push_back(i);
            y.push_back(*iter);
            val.push_back(rand() % (max_value - min_value) + min_value);
        }
    }
}

template <typename I, typename T>
inline bool compare(const std::tuple<I, I, T, I> &lhs, const std::tuple<I, I, T, I> &rhs) {
    I a = std::get<0>(lhs);
    I b = std::get<0>(rhs);
    I c = std::get<1>(lhs);
    I d = std::get<1>(rhs);
    if (a == b)
        return c < d;
    else
        return a < b;
}

template <typename I>
inline bool compare(const std::tuple<I, I, I> &lhs, const std::tuple<I, I, I> &rhs) {
    I a = std::get<0>(lhs);
    I b = std::get<0>(rhs);
    I c = std::get<1>(lhs);
    I d = std::get<1>(rhs);
    if (a == b)
        return c < d;
    else
        return a < b;
}

template <typename I, typename T>
inline void customSort(I *row_indices, I *col_indices, T *values, I nnz) {
    I nvals = nnz;
    std::vector<std::tuple<I, I, T, I>> my_tuple;

    for (I i = 0; i < nvals; ++i)
        my_tuple.push_back(std::make_tuple(row_indices[i], col_indices[i], values[i], i));

    std::sort(my_tuple.begin(), my_tuple.end(), compare<I, T>);

    for (I i = 0; i < nvals; ++i) {
        row_indices[i] = std::get<0>(my_tuple[i]);
        col_indices[i] = std::get<1>(my_tuple[i]);
        values[i] = std::get<2>(my_tuple[i]);
    }
}

template <typename I, typename T>
inline void customSort(std::vector<I> *row_indices, std::vector<I> *col_indices, std::vector<T> *values) {
    customSort(row_indices->data(), col_indices->data(), values->data(), (I) row_indices->size());
}

template <typename I, typename T>
inline void coo2csr(I *csrRowPtr, I *csrColInd, T *csrVal, I *row_indices,
                    I *col_indices, T *values, I nrows, I ncols, I nnz) {

    I temp, row, col, dest, cumsum = 0;

    std::vector<I> row_indices_t(row_indices, row_indices + nnz);
    std::vector<I> col_indices_t(col_indices, col_indices + nnz);
    std::vector<T> values_t(values, values + nnz);

    customSort<I, T>(row_indices_t.data(), col_indices_t.data(), values_t.data(), nnz);

    // Set all rowPtr to 0
    for (I i = 0; i <= nrows; i++)
        csrRowPtr[i] = 0;

    // Go through all elements to see how many fall in each row
    for (I i = 0; i < nnz; i++) {
        row = row_indices_t[i];
        if (row >= nrows)
            std::cout << "Error: Index out of bounds!\n";
        csrRowPtr[row]++;
    }

    // Cumulative sum to obtain rowPtr
    for (I i = 0; i < nrows; i++) {
        temp = csrRowPtr[i];
        csrRowPtr[i] = cumsum;
        cumsum += temp;
    }
    csrRowPtr[nrows] = nnz;

    // Store colInd and val
    for (I i = 0; i < nnz; i++) {
        row = row_indices_t[i];
        dest = csrRowPtr[row];
        col = col_indices_t[i];
        if (col >= ncols)
            std::cout << "Error: Index out of bounds!\n";
        csrColInd[dest] = col;
        csrVal[dest] = values_t[i];
        csrRowPtr[row]++;
    }
    cumsum = 0;

    // Undo damage done to rowPtr
    for (I i = 0; i < nrows; i++) {
        temp = csrRowPtr[i];
        csrRowPtr[i] = cumsum;
        cumsum = temp;
    }
    temp = csrRowPtr[nrows];
    csrRowPtr[nrows] = cumsum;
    cumsum = temp;
}

template <typename T>
inline void print_graph(T *ptr, T *idx, int N, int max_N = 20, int max_E = 20) {
    std::cout << "-) degree: " << ptr[0] << std::endl;
    for (int v = 1; v < std::min((int)N + 1, max_N); v++) {
        std::cout << v - 1 << ") degree: " << ptr[v] - ptr[v - 1] << ", edges: ";
        for (int e = 0; e < ptr[v] - ptr[v - 1]; e++) {
            if (e < max_E) {
                std::cout << idx[ptr[v - 1] + e] << ", ";
            }
        }
        std::cout << std::endl;
    }
}

inline std::string random_string(const int len) {
    std::string tmp_s;
    static const char ALPHANUM[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    
    for (int i = 0; i < len; i++) {
        tmp_s += ALPHANUM[rand() % (sizeof(ALPHANUM) - 1)];
    }
    return tmp_s;
}

inline void print_string(std::string s, int max_print = 20) {
    if (s.length() > max_print) {
        std::cout << s.substr(0, max_print / 2) << "..." << s.substr(s.length() - max_print / 2 - 1, s.length()) << std::endl;
    } else {
        std::cout << s << std::endl;
    }
}

template<typename T>
inline int check_array_equality(T *x, T *y, int n, float tol = 0.0000001f, bool debug = false, int max_print = 20) {
	int num_errors = 0;
	for (int i = 0; i < n; i++) {
		T diff = (T)((x[i] > y[i]) ? (x[i] - y[i]) : (y[i] - x[i]));
		if ((float) diff > tol) {
			num_errors++;
			if (debug && num_errors < max_print) {
				std::cout << i << ") X: " << x[i] << ", Y: " << y[i] << ", diff: " << diff << std::endl;
			}
		}
	}
	return num_errors;
}

template<typename T>
inline int check_array_equality(std::vector<T> &x, std::vector<T> &y, float tol = 0.0000001f, bool debug = false, int max_print = 20) {
    return check_array_equality(x.data(), y.data(), x.size(), tol, debug, max_print);
}

template<typename T>
inline bool check_equality(T x, T y, float tol = 0.0000001f, bool debug = false) {
	bool equal = true;

	float diff = std::abs(x - y);
	if (diff > tol) {
		equal = false;
		if (debug) {
			std::cout << "x: " << x << ", y: " << y << ", diff: " << diff << std::endl;
		}
	}
	return equal;
}

///////////////////////////////
///////////////////////////////

// Utility functions adapted from Graphblast, used to read MTX files;

template<typename I, typename T>
inline void readTuples(std::vector<I> *row_indices, std::vector<I> *col_indices,
		std::vector<T> *values, I nvals, FILE *f, bool read_values = true,
		bool zero_indexed_file = false) {
	I row_ind, col_ind;
	int_type row_ind_i, col_ind_i;
	double value;
	double acc = 0.0;
	std::vector<double> tmp_value_store;
	// Currently checks if there are fewer rows than promised
	// Could add check for edges in diagonal of adjacency matrix
	for (I i = 0; i < nvals; i++) {
		if (fscanf(f, "%u", &row_ind_i) == EOF) {
//			std::cout << "Error: Not enough rows in mtx file!\n";
//			return;
			row_ind = i;
			col_ind = 1;
			value = (T) 0.0;
		} else {
			fscanf(f, "%u", &col_ind_i);
			if (read_values) {
				fscanf(f, "%lf", &value);
			} else {
				value = 1;
			}
			row_ind = (I) row_ind_i;
			col_ind = (I) col_ind_i;

			// Convert 1-based indexing MTX to 0-based indexing C++
			if (!zero_indexed_file) {
				row_ind--;
				col_ind--;
			}
			row_indices->push_back(row_ind);
			col_indices->push_back(col_ind);
			tmp_value_store.push_back(value);
			double d_value = double(value);
			acc += d_value * d_value;
		}
	}

	for (I i = 0; i < nvals; ++i) {
		double d_value = tmp_value_store[i];
		double d_res = d_value / acc;
		values->push_back(T(std::sqrt(d_res)));
	}

}

template<typename I, typename T>
inline void undirect(std::vector<I> *row_indices, std::vector<I> *col_indices,
		std::vector<T> *values, I *nvals,
		bool skip_values = false) {
	for (I i = 0; i < *nvals; i++) {
		if ((*col_indices)[i] != (*row_indices)[i]) {
			row_indices->push_back((*col_indices)[i]);
			col_indices->push_back((*row_indices)[i]);
			if (!skip_values) {
				values->push_back((*values)[i]);
			}
		}
	}
	*nvals = row_indices->size();
}

/*!
 * Remove self-loops, duplicates and make graph undirected if option is set
 */
template<typename I, typename T>
inline void removeSelfloop(std::vector<I> *row_indices,
		std::vector<I> *col_indices, std::vector<T> *values, I *nvals) {
	// Sort
	customSort<I, T>(row_indices, col_indices, values);

	I curr = (*col_indices)[0];
	I last;
	I curr_row = (*row_indices)[0];
	I last_row;

	// Detect self-loops and duplicates
	for (I i = 0; i < *nvals; i++) {
		last = curr;
		last_row = curr_row;
		curr = (*col_indices)[i];
		curr_row = (*row_indices)[i];

		// Self-loops
		if (curr_row == curr)
			(*col_indices)[i] = -1;

		// Duplicates
		if (i > 0 && curr == last && curr_row == last_row)
			(*col_indices)[i] = -1;
	}

	I shift = 0;

	// Remove self-loops and duplicates marked -1.
	I back = 0;
	for (I i = 0; i + shift < *nvals; i++) {
		if ((*col_indices)[i] == -1) {
			for (; back <= *nvals; shift++) {
				back = i + shift;
				if ((*col_indices)[back] != -1) {
					(*col_indices)[i] = (*col_indices)[back];
					(*row_indices)[i] = (*row_indices)[back];
					(*col_indices)[back] = -1;
					break;
				}
			}
		}
	}

	*nvals = *nvals - shift;
	row_indices->resize(*nvals);
	col_indices->resize(*nvals);
	values->resize(*nvals);
}

template<typename I, typename T>
inline int readMtx(const char *fname, std::vector<I> *row_indices,
		std::vector<I> *col_indices, std::vector<T> *values, I* num_rows,
		I* num_cols, I* num_nnz, int directed = true, bool read_values = false,
		bool debug = false, bool zero_indexed_file = false, bool sort_tuples =
				false) {
	int ret_code;
	MM_typecode matcode;
	FILE *f;

	I nrows = 0;
	I ncols = 0;
	I nvals = 0;
	bool mtxinfo = false;

	if ((f = fopen(fname, "r")) == NULL) {
		std::cerr << "File " << fname << " not found" << std::endl;
		std::cerr.flush();
		exit(1);
	}

	// Read MTX banner
	if (mm_read_banner(f, &matcode) != 0) {
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	// Read MTX Size
	if ((ret_code = mm_read_mtx_crd_size(f, &nrows, &ncols, &nvals)) != 0)
		exit(1);
	readTuples<I, T>(row_indices, col_indices, values, nvals, f, read_values,
			zero_indexed_file);

	bool is_undirected = mm_is_symmetric(matcode) || directed == 2;
	is_undirected = (directed == 1) ? false : is_undirected;
	if (is_undirected) {
		undirect(row_indices, col_indices, values, &nvals);
	}
	if (sort_tuples)
		customSort<I, T>(row_indices, col_indices, values);

	if (mtxinfo)
		mm_write_banner(stdout, matcode);
	if (mtxinfo)
		mm_write_mtx_crd_size(stdout, nrows, ncols, nvals);

	*num_rows = nrows;
	*num_cols = ncols;
	*num_nnz = nvals;

	return ret_code;
}