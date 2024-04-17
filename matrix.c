#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/ 

/*
 * Generates a random double between `low` and `high`.
 */
double rand_double(double low, double high) {
    double range = (high - low); 
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/*
 * Generates a random matrix with `seed`.
 */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocate space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. Remember to set all fields of the matrix struct.
 * `parent` should be set to NULL to indicate that this matrix is not a slice.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    
    // Check validity of dimensions, allocate data memory
    if (rows < 1 || cols < 1) {
	PyErr_SetString(PyExc_ValueError, "Invalid inputs");
        return -1;
    }

    double *d_doubles = (double*) calloc(rows * cols, sizeof(double));
    if (d_doubles == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
	    return -1;
    }
    // double **d = (double**) calloc(rows, sizeof(double*));
    // if (d == NULL) {
    //     PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
	//     return -1;
    // }

    double **d = (double**) calloc(rows, sizeof(double*));
    if (d == NULL) {
	free(d_doubles);
        PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
	    return -1;
    }
    for (int i = 0; i < rows; i++) {
        d[i] = d_doubles + (i * cols);
    }

    // for (int i = 0; i < rows; i++) {
    //     d[i] = (double*) calloc(cols, sizeof(double));
    //     if (d[i] == NULL) {
    //         for (int j = 0; j < i; j++) {
    //             free(d[j]);
    //         }
	//         PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
    //         return -1;
    //     }
    // }

    // Allocate struct memory
    matrix *m = (matrix*) malloc(sizeof(matrix));
    if (m == NULL) {
        // for (int j = 0; j < rows; j++) {
  	    //     free(d[j]);
        // }
	free(d_doubles);
        free(d);
	PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
        return -1;
    }
    m->rows = rows;
    m->cols = cols;

    m->data = d;
    m->double_data = d_doubles;

    if (rows == 1 || cols == 1) {
        m->is_1d = 1;
    } else {
        m->is_1d = 0;
    }
    m->ref_cnt = 1;
    m->parent = NULL;

    // Set passed-in pointer to allocated memory
    *mat = m;

    return 0;
    
}

/*
 * Allocate space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * This is equivalent to setting the new matrix to be
 * from[row_offset:row_offset + rows, col_offset:col_offset + cols]
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                        int rows, int cols) {

    // Check validity of dimensions, allocate data memory
    if (rows < 1 || cols < 1) {
        PyErr_SetString(PyExc_ValueError, "Invalid inputs");
	    return -1;
    }
    if (from->rows < rows || from->cols < cols) {
        PyErr_SetString(PyExc_IndexError, "Index out of range");
	    return -1;
    }

    double **d = (double**) calloc(rows, sizeof(double*));
    if (d == NULL) {
	    PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
        return -1;
    }

    if (cols != 1) {
	# pragma omp parallel for 
        for (int i = 0; i < rows; i++) {
            d[i] = &((from->data)[i + row_offset][col_offset]); // what is the resultant type of d[i] and width of memory?
        }
    } else {
	# pragma omp parallel for 
	for (int i = 0; i < rows; i++) {
            d[i] = &((from->data)[i + row_offset][col_offset]); // what is the resultant type of d[i] and width of memory?
        }
    }

    // Allocate struct memory
    matrix *m = (matrix*) malloc(sizeof(matrix));
    if (m == NULL) {
	free(d);
        PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
	return -1;
    }

    m->rows = rows;
    m->cols = cols;
    m->data = d;
    m->double_data = from->double_data; // hopefully we shouldn't actually have to use this...

    if (rows == 1 || cols == 1) {
        m->is_1d = 1;
    } else {
        m->is_1d = 0;
    }
	
    matrix *parent_ptr = from;
    while (parent_ptr->parent != NULL) {
        parent_ptr = parent_ptr->parent;
    }
    parent_ptr->ref_cnt += 1; // FLAG AS DANGEROUS
    m->parent = parent_ptr;

    m->ref_cnt = 1; // does not matter what this is, since m is not, nor will ever be, the root

    // Set passed-in pointer to allocated memory
    *mat = m;

    return 0;


	// int rtval = allocate_matrix(mat, rows, cols);   
	// if (rtval != 0) {
	// 	return -1;
	// }
	// for (int i = row_offset; i < row_offset + rows; i++) {
	// 	for (int j = col_offset; j < col_offset + cols; j++) {
	// 		((*mat)->data)[i - row_offset][j - col_offset] = (from->data)[i][j];
	// 	}
	// }
	// (*mat)->parent = from;
	// from->ref_cnt += 1;
	// return 0;
}

/*
 * This function will be called automatically by Python when a numc matrix loses all of its
 * reference pointers.
 * You need to make sure that you only free `mat->data` if no other existing matrices are also
 * referring this data array.
 * See the spec for more information.
 */
void deallocate_matrix(matrix *mat) {
    if (mat == NULL) {
        return;
    }
    // if we are not at the root
    if (mat->parent != NULL) {
        deallocate_matrix(mat->parent);
	free(mat->data);
        free(mat);
    } else if (mat->parent == NULL) {
        if (mat->ref_cnt == 1) {
            // for (int i = 0; i < mat->rows; i++) { 
            //     free((mat->data)[i]);
            // }
            free(mat->double_data);
            free(mat->data);
            free(mat);
        } else {
            mat->ref_cnt -= 1;
        }
    } else {
        PyErr_SetString(PyExc_RuntimeError, "Shouldn't get to this case of deallocate_matrix");
    }

    // // don't do anything if it has a parent?
    // if (mat->parent != NULL) {
    //     // deallocate matrix ref work
    //     // ensure that parent is alwayys real matrix not other ref
    //     // check if from matrix has a parent --> use from
    //     // only deallocate for the parent if all parents and refs are none
    //     return;
    // }
    // // if count of the from matrix goes to zero then deallcoate the parents
    // if (mat->ref_cnt != 0) {
    //     return;
    // }
    // // deallocate all the {} data in the matrix
    // for (int i = 0; i < mat->rows; i++) {
    //     free((mat->data)[i]);
    // }
    // free(mat->data);
    // // deallocate the matrix struct itself
    // free(mat);
}

/*
 * Return the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    return (mat->data)[row][col];
}

/*
 * Set the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    (mat->data)[row][col] = val;
}

/*
 * Set all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    if (mat == NULL) {
        return;
    }
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; i < mat->cols; j++) {
            (mat->data)[i][j] = val;
        }
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    if (mat1->rows != mat2->rows) {
        PyErr_SetString(PyExc_ValueError, "Invalid input");
	    return -1;
    }
    if (mat1->cols != mat2->cols) {
	    PyErr_SetString(PyExc_ValueError, "Invalid input");
        return -1;
    }

    if (mat1->rows * mat1->cols > 50000) {
        int N = mat1->rows * mat1->cols;

        # pragma omp parallel for 
        for (int i = 0; i < ((N / 4) * 4) ; i += 4) {
            __m256d vector1 = _mm256_loadu_pd((__m256d*) (mat1->double_data + i));
            __m256d vector2 = _mm256_loadu_pd((__m256d*) (mat2->double_data + i));
            __m256d sum = _mm256_add_pd(vector1, vector2);
            // printf("sum vector was %d \n", sum); // putting this print statement caused an infinite loop of printing the same number...
            _mm256_storeu_pd(result->double_data + i, sum);
            if (i == N/4) {
                printf("number of threads used was %d \n", omp_get_num_threads());
            }
        }
        
        for (int i = (N / 4) * 4; i < N; i++) {
            *((result->double_data) + i) = *((mat1->double_data) + i) + *((mat2->double_data) + i); 
        }
    } else {
        // # pragma omp parallel for 
        for (int i = 0; i < mat1->rows * mat1->cols; i++) {
            *((result->double_data) + i) = *((mat1->double_data) + i) + *((mat2->double_data) + i);
        }

        // for (int i = 0; i < mat1->rows; i++) {
        //     for (int j = 0; j < mat2->cols; j++) {
        //         (result->data)[i][j] = (mat1->data)[i][j] + (mat2->data)[i][j];
        //     }
        // }   
    }
    
    return 0;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    if (mat1->rows != mat2->rows) {
        PyErr_SetString(PyExc_ValueError, "Invalid input");
	    return -1;
    }
    if (mat1->cols != mat2->cols) { 
	    PyErr_SetString(PyExc_ValueError, "Invalid input");
        return -1; 
    }
     
    if (mat1->rows * mat1->cols > 50000) {
        int N = mat1->rows * mat1->cols;

        # pragma omp parallel for 
        for (int i = 0; i < ((N / 4) * 4) ; i += 4) {
            __m256d vector1 = _mm256_loadu_pd((__m256d*) (mat1->double_data + i));
            __m256d vector2 = _mm256_loadu_pd((__m256d*) (mat2->double_data + i));
            __m256d diff = _mm256_sub_pd(vector1, vector2);
            // printf("sum vector was %d \n", sum); // putting this print statement caused an infinite loop of printing the same number...
            _mm256_storeu_pd(result->double_data + i, diff);
            if (i == N/4) {
                printf("number of threads used was %d \n", omp_get_num_threads());
            }
        }

        for (int i = (N / 4) * 4; i < N; i++) {
            *((result->double_data) + i) = *((mat1->double_data) + i) - *((mat2->double_data) + i);
        }
    } else {
        // # pragma omp parallel for 
        for (int i = 0; i < mat1->rows * mat1->cols; i++) {
            *((result->double_data) + i) = *((mat1->double_data) + i) - *((mat2->double_data) + i);
        }
    }

    
   // for (int i = 0; i < mat1->rows; i++) {
   //     for (int j = 0; j < mat2->cols; j++) {
   //         (result->data)[i][j] = (mat1->data)[i][j] - (mat2->data)[i][j];
   //     }
   // }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    if (mat1->cols != mat2->rows) {
	    PyErr_SetString(PyExc_ValueError, "Invalid input");
        return -1;
    }

    double* trans = NULL;
    // trans = (double*) calloc(mat2->rows * mat2->cols, sizeof(double));
    // if (trans == NULL) {
	//     PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
    //     return -1;
    // }
    // for (int i = 0; i < mat2->rows * mat2->cols; i++) {
    //     *(trans + (i % mat2->cols)* mat2->rows + (i / mat2->cols)) = mat2->double_data[i];
    // }

    trans = (double*) calloc(mat2->rows * mat2->cols, sizeof(double));
    # pragma omp parallel for 
    for (int i = 0; i < mat2->rows * mat2->cols; i++) {
        *(trans + (i % mat2->cols)* mat2->rows + (i / mat2->cols)) = mat2->double_data[i];
    }
        // for (int i = 0; i < mat2->rows * mat2->cols; i++) {
        //     *(trans + i) = mat2->data[i % mat2->rows][i / mat2->rows];
        // }
    
    // matrix *mat2T = NULL;
    // allocate_matrix(&mat2T, mat2->cols, mat2->rows);
    // for (int i = 0; i < mat2->cols; i++) {
    //     for (int j = 0; j < mat2->rows; j++) {
    //         mat2T->data[i][j] = mat2->data[j][i];
    //     }
    // }

    // for (int i = 0; i < mat2->cols; i++) {
    //     for (int j = 0; j < mat2->rows; j++) {
    //         *(d + j*mat2->cols + i) = mat2->data[j][i];
    //     }
    // }

    if (mat1->rows * mat1->cols * mat2->cols > 20000) {
        # pragma omp parallel for 
        for (int r1 = 0; r1 < mat1->rows; r1++) {
            for (int r2 = 0; r2 < mat2->cols; r2++) {
                
                double *m1 = (mat1->data)[r1];
                double *m2 = trans + r2*mat2->rows;
                
                // for (int j = 0; j < mat1->cols; j++) {
                //     sum += m1[j] * m2[j];
                // }
                // (result->data)[r1][r2] = sum;
 			
                __m256d suma = _mm256_set1_pd(0); 
		__m256d sumb = _mm256_set1_pd(0);
		__m256d sumc = _mm256_set1_pd(0);
		__m256d sumd = _mm256_set1_pd(0);

                for (int i = 0; i < ((mat1->cols / 16) * 16) ; i += 16) {
                    __m256d vector1 = _mm256_loadu_pd((__m256d*) (m1 + i));
                    __m256d vector2 = _mm256_loadu_pd((__m256d*) (m2 + i));
		    suma = _mm256_fmadd_pd(vector1, vector2, suma);
		     		
		    vector1 = _mm256_loadu_pd((__m256d*) (m1 + i + 4));
                    vector2 = _mm256_loadu_pd((__m256d*) (m2 + i + 4));
                    sumb = _mm256_fmadd_pd (vector1, vector2, sumb);
		
		    vector1 = _mm256_loadu_pd((__m256d*) (m1 + i + 8));
                    vector2 = _mm256_loadu_pd((__m256d*) (m2 + i + 8));
                    sumc = _mm256_fmadd_pd(vector1, vector2, sumc);

		    vector1 = _mm256_loadu_pd((__m256d*) (m1 + i + 12));
                    vector2 = _mm256_loadu_pd((__m256d*) (m2 + i + 12));
                    sumd = _mm256_fmadd_pd(vector1, vector2, sumd);
		}

		for (int i = ((mat1->cols / 16) * 16); i < ((mat1->cols / 4) * 4); i += 4) {
		    __m256d vector1 = _mm256_loadu_pd((__m256d*) (m1 + i));
                    __m256d vector2 = _mm256_loadu_pd((__m256d*) (m2 + i));
                    suma = _mm256_fmadd_pd(vector1, vector2, suma);
		}
		
		suma = _mm256_add_pd(sumb, suma);
		suma = _mm256_add_pd(sumc, suma);
		suma = _mm256_add_pd(sumd, suma);
		
		double resulta[4];
                _mm256_storeu_pd((__m256d*) resulta, suma);
                double final_sum = resulta[0] + resulta[1] + resulta[2] + resulta[3];

                for (int i = ((mat1->cols / 4) * 4); i < mat1->cols; i++) {
                    final_sum += *(m1 + i) * *(m2 + i); 
                }

                (result->data)[r1][r2] = final_sum;
            }
        }
    } else {
        for (int r1 = 0; r1 < mat1->rows; r1++) {
            for (int r2 = 0; r2 < mat2->cols; r2++) {
                double sum = 0;
                double *m1 = (mat1->data)[r1];
                double *m2 = trans + r2*mat2->rows;
                for (int j = 0; j < mat1->cols; j++) {
                    sum += m1[j] * m2[j];
                }
                (result->data)[r1][r2] = sum;
            }
        }
    }

    // deallocate_matrix(mat2T);

    // for (int i = 0; i < mat1->rows; i++) {
    //     for (int j = 0; j < mat2->cols; j++) {
    //         double sum = 0;
    //         for (int k = 0; k < mat1->cols; k++) {
    //             sum += (mat1->data)[i][k] * (mat2->data)[k][j];
    //         }
    //         (result->data)[i][j] = sum;
    //     }
    // }

    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    if (pow == 0) {
        # pragma omp parallel for
        for (int i = 0; i < mat->rows; i++) {
	    (result->data)[i][i] = 1;
	} 
	return 0;
    }

    if (pow == 1) {
	# pragma omp parallel for 
	for (int i = 0; i < mat->rows; i++) {
            for (int j = 0; j < mat->cols; j++) {
                (result->data)[i][j] = (mat->data)[i][j];
            }
        }
        return 0;
     }

     matrix *square = NULL;
     int h = allocate_matrix(&square, mat->rows, mat->cols);
     if (h != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
        return -1;
     }

     mul_matrix(square, mat, mat);

     if (pow - 2*(pow/2) == 0) {
	pow_matrix(result, square, pow/2);
     } else {
        matrix *intermediate_result = NULL;
        int h = allocate_matrix(&intermediate_result, mat->rows, mat->cols);
        if (h != 0) {
	    PyErr_SetString(PyExc_RuntimeError, "Ran out of memory");
            return -1;
        }
        pow_matrix(intermediate_result, square, pow/2);
        mul_matrix(result, intermediate_result, mat);
        deallocate_matrix(intermediate_result);
     }

     deallocate_matrix(square);

     return 0;

}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    
    int N = mat->rows * mat->cols;
    __m256d zeros = _mm256_set1_pd(0); 

    # pragma omp parallel for 
    for (int i = 0; i < ((N / 4) * 4) ; i += 4) {
        __m256d vector = _mm256_loadu_pd((__m256d*) (mat->double_data + i));
        __m256d neg_vector = _mm256_sub_pd(zeros, vector);
        _mm256_storeu_pd(result->double_data + i, neg_vector);
    }
    
    for (int i = (N / 4) * 4; i < N; i++) {
        *((result->double_data) + i) = 0.0 - *((mat->double_data) + i); 
    }

    // for (int i = 0; i < mat->rows; i++) {
    //     for (int j = 0; j < mat->cols; j++) {
    //         (result->data)[i][j] = -(mat->data)[i][j];
    //     }
    // }

    return 0;
}
/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */

static double double_mask = nan("0xfffffffffffff");

int abs_matrix(matrix *result, matrix *mat) {

    int N = mat->rows * mat->cols;

    __m256d mask = _mm256_set1_pd(double_mask); 

    # pragma omp parallel for 
    for (int i = 0; i < ((N / 4) * 4) ; i += 4) {
        __m256d vector = _mm256_loadu_pd((__m256d*) (mat->double_data + i));
        __m256d abs_vector = _mm256_and_pd(mask, vector);
        _mm256_storeu_pd(result->double_data + i, abs_vector);
    }
    
    for (int i = (N / 4) * 4; i < N; i++) {
        if (*((mat->double_data) + i) < 0) {
		    *((result->double_data) + i) = 0.0 - *((mat->double_data) + i); 
	    } else {
            *(result->double_data + i) = *((mat->double_data) + i);
	    }
    }

    // for (int i = 0; i < mat->rows; i++) {
    //     for (int j = 0; j < mat->cols; j++) {
    //         (result->data)[i][j] = -(mat->data)[i][j];
    //     }
    // }

    // for (int i = 0; i < mat->rows; i++) {
    //     for (int j = 0; j < mat->cols; j++) {
	//     double value = (mat->data)[i][j];
	//     if (value < 0) {
	// 	    (result->data)[i][j] = 0.0 - value;
	//     } else {
    //         (result->data)[i][j] = value;
	//     }
	// }
    // }
    return 0;
}
