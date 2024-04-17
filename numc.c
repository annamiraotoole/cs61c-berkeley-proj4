#include "numc.h"
#include <structmember.h>

PyTypeObject Matrix61cType;

/* Helper functions for initalization of matrices and vectors */

/*
 * Return a tuple given rows and cols
 */
PyObject *get_shape(int rows, int cols) {
  if (rows == 1 || cols == 1) {
    return PyTuple_Pack(1, PyLong_FromLong(rows * cols));
  } else {
    return PyTuple_Pack(2, PyLong_FromLong(rows), PyLong_FromLong(cols));
  }
}
/*
 * Matrix(rows, cols, low, high). Fill a matrix random double values
 */
int init_rand(PyObject *self, int rows, int cols, unsigned int seed, double low,
              double high) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    rand_matrix(new_mat, seed, low, high);
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(rows, cols, val). Fill a matrix of dimension rows * cols with val
 */
int init_fill(PyObject *self, int rows, int cols, double val) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed)
        return alloc_failed;
    else {
        fill_matrix(new_mat, val);
        ((Matrix61c *)self)->mat = new_mat;
        ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    }
    return 0;
}

/*
 * Matrix(rows, cols, 1d_list). Fill a matrix with dimension rows * cols with 1d_list values
 */
int init_1d(PyObject *self, int rows, int cols, PyObject *lst) {
    if (rows * cols != PyList_Size(lst)) {
        PyErr_SetString(PyExc_ValueError, "Incorrect number of elements in list");
        return -1;
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j, PyFloat_AsDouble(PyList_GetItem(lst, count)));
            count++;
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(2d_list). Fill a matrix with dimension len(2d_list) * len(2d_list[0])
 */
int init_2d(PyObject *self, PyObject *lst) {
    int rows = PyList_Size(lst);
    if (rows == 0) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot initialize numc.Matrix with an empty list");
        return -1;
    }
    int cols;
    if (!PyList_Check(PyList_GetItem(lst, 0))) {
        PyErr_SetString(PyExc_ValueError, "List values not valid");
        return -1;
    } else {
        cols = PyList_Size(PyList_GetItem(lst, 0));
    }
    for (int i = 0; i < rows; i++) {
        if (!PyList_Check(PyList_GetItem(lst, i)) ||
                PyList_Size(PyList_GetItem(lst, i)) != cols) {
            PyErr_SetString(PyExc_ValueError, "List values not valid");
            return -1;
        }
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j,
                PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(lst, i), j)));
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * This deallocation function is called when reference count is 0
 */
void Matrix61c_dealloc(Matrix61c *self) {
    deallocate_matrix(self->mat);
    Py_TYPE(self)->tp_free(self);
}

/* For immutable types all initializations should take place in tp_new */
PyObject *Matrix61c_new(PyTypeObject *type, PyObject *args,
                        PyObject *kwds) {
    /* size of allocated memory is tp_basicsize + nitems*tp_itemsize*/
    Matrix61c *self = (Matrix61c *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

/*
 * This matrix61c type is mutable, so needs init function. Return 0 on success otherwise -1
 */
int Matrix61c_init(PyObject *self, PyObject *args, PyObject *kwds) {
    /* Generate random matrices */
    if (kwds != NULL) {
        PyObject *rand = PyDict_GetItemString(kwds, "rand");
        if (!rand) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (!PyBool_Check(rand)) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (rand != Py_True) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        PyObject *low = PyDict_GetItemString(kwds, "low");
        PyObject *high = PyDict_GetItemString(kwds, "high");
        PyObject *seed = PyDict_GetItemString(kwds, "seed");
        double double_low = 0;
        double double_high = 1;
        unsigned int unsigned_seed = 0;

        if (low) {
            if (PyFloat_Check(low)) {
                double_low = PyFloat_AsDouble(low);
            } else if (PyLong_Check(low)) {
                double_low = PyLong_AsLong(low);
            }
        }

        if (high) {
            if (PyFloat_Check(high)) {
                double_high = PyFloat_AsDouble(high);
            } else if (PyLong_Check(high)) {
                double_high = PyLong_AsLong(high);
            }
        }

        if (double_low >= double_high) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        // Set seed if argument exists
        if (seed) {
            if (PyLong_Check(seed)) {
                unsigned_seed = PyLong_AsUnsignedLong(seed);
            }
        }

        PyObject *rows = NULL;
        PyObject *cols = NULL;
        if (PyArg_UnpackTuple(args, "args", 2, 2, &rows, &cols)) {
            if (rows && cols && PyLong_Check(rows) && PyLong_Check(cols)) {
                return init_rand(self, PyLong_AsLong(rows), PyLong_AsLong(cols), unsigned_seed, double_low,
                                 double_high);
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    }
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 3, &arg1, &arg2, &arg3)) {
        /* arguments are (rows, cols, val) */
        if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && (PyLong_Check(arg3)
                || PyFloat_Check(arg3))) {
            if (PyLong_Check(arg3)) {
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyLong_AsLong(arg3));
            } else
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyFloat_AsDouble(arg3));
        } else if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && PyList_Check(arg3)) {
            /* Matrix(rows, cols, 1D list) */
            return init_1d(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), arg3);
        } else if (arg1 && PyList_Check(arg1) && arg2 == NULL && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_2d(self, arg1);
        } else if (arg1 && arg2 && PyLong_Check(arg1) && PyLong_Check(arg2) && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), 0);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return -1;
    }
}

/*
 * List of lists representations for matrices
 */
PyObject *Matrix61c_to_list(Matrix61c *self) {
    int rows = self->mat->rows;
    int cols = self->mat->cols;
    PyObject *py_lst = NULL;
    if (self->mat->is_1d) {  // If 1D matrix, print as a single list
        py_lst = PyList_New(rows * cols);
        int count = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(py_lst, count, PyFloat_FromDouble(get(self->mat, i, j)));
                count++;
            }
        }
    } else {  // if 2D, print as nested list
        py_lst = PyList_New(rows);
        for (int i = 0; i < rows; i++) {
            PyList_SetItem(py_lst, i, PyList_New(cols));
            PyObject *curr_row = PyList_GetItem(py_lst, i);
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(curr_row, j, PyFloat_FromDouble(get(self->mat, i, j)));
            }
        }
    }
    return py_lst;
}

PyObject *Matrix61c_class_to_list(Matrix61c *self, PyObject *args) {
    PyObject *mat = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 1, &mat)) {
        if (!PyObject_TypeCheck(mat, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        Matrix61c* mat61c = (Matrix61c*)mat;
        return Matrix61c_to_list(mat61c);
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
}

/*
 * Add class methods
 */
PyMethodDef Matrix61c_class_methods[] = {
    {"to_list", (PyCFunction)Matrix61c_class_to_list, METH_VARARGS, "Returns a list representation of numc.Matrix"},
    {NULL, NULL, 0, NULL}
};

/*
 * Matrix61c string representation. For printing purposes.
 */
PyObject *Matrix61c_repr(PyObject *self) {
    PyObject *py_lst = Matrix61c_to_list((Matrix61c *)self);
    return PyObject_Repr(py_lst);
}

/* NUMBER METHODS */

Matrix61c *number_method_unpack(PyObject* args) {
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    Matrix61c* mat61c_b = (Matrix61c*) args;
    return mat61c_b;
}

int check_matrix61c_dims(Matrix61c* a, Matrix61c* b) {
    if (b->mat == NULL) {
        return -1;
    }
    if (a->mat == NULL) {
        return -1;
    }
    if (b->mat->rows != a->mat->rows) {
        PyErr_SetString(PyExc_ValueError, "Incorrect matrix dimensions");
        return -1;
    }
    if (b->mat->cols != a->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Incorrect matrix dimensions");
        return -1;
    }
    return 0;
}


/*
 * Add the second numc.Matrix (Matrix61c) object to the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_add(Matrix61c* self, PyObject* args) {
    Matrix61c* b = number_method_unpack(args);
    if (b == NULL) {
        return NULL;
    }
    // check dimensions match for a and b
    if (check_matrix61c_dims(self, b)) {
        return NULL;
    }

    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    matrix *new_mat = NULL;
    allocate_matrix(&new_mat, b->mat->rows, b->mat->cols);
    result->mat = new_mat;
    result->shape = get_shape(b->mat->rows, b->mat->cols);

    int retval = add_matrix(result->mat, self->mat, b->mat);
    if (retval != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix addition failed");
        return NULL;
    }

    return (PyObject *) result;
}

/*
 * Substract the second numc.Matrix (Matrix61c) object from the first one. The first operand is 
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_sub(Matrix61c* self, PyObject* args) {
    Matrix61c* b = number_method_unpack(args);
    if (b == NULL) {
        return NULL;
    }
    // check dimensions match for a and b
    if (check_matrix61c_dims(self, b)) {
        return NULL;
    }

    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    matrix *new_mat = NULL;
    allocate_matrix(&new_mat, b->mat->rows, b->mat->cols);
    result->mat = new_mat;
    result->shape = get_shape(b->mat->rows, b->mat->cols);

    int retval = sub_matrix(result->mat, self->mat, b->mat);
    if (retval != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix subtraction failed");
        return NULL;
    }

    return (PyObject *) result;
}

/*
 * NOT element-wise multiplication. The first operand is self, and the second operand
 * can be obtained by casting `args`.
 */
PyObject *Matrix61c_multiply(Matrix61c* self, PyObject *args) {
    Matrix61c* b = number_method_unpack(args);
    if (b == NULL) {
        return NULL;
    }
    // check dimensions match for a and b
    if (self->mat->cols != b->mat->rows) {
        PyErr_SetString(PyExc_ValueError, "Incorrect matrix dimensions");
        return NULL;
    }

    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    matrix *new_mat = NULL;
    allocate_matrix(&new_mat, self->mat->rows, b->mat->cols);
    result->mat = new_mat;
    result->shape = get_shape(self->mat->rows, b->mat->cols);

    int retval = mul_matrix(result->mat, self->mat, b->mat);
    if (retval != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix multiplication failed");
        return NULL;
    }

    return (PyObject *) result;
}

/*
 * Negates the given numc.Matrix.
 */
PyObject *Matrix61c_neg(Matrix61c* self) {
    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    matrix *new_mat = NULL;
    allocate_matrix(&new_mat, self->mat->rows, self->mat->cols);
    result->mat = new_mat;
    result->shape = get_shape(self->mat->rows, self->mat->cols);

    int retval = neg_matrix(result->mat, self->mat);
    if (retval != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix negation failed");
        return NULL;
    }

    return (PyObject *) result;
}

/*
 * Take the element-wise absolute value of this numc.Matrix.
 */
PyObject *Matrix61c_abs(Matrix61c *self) {
    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    matrix *new_mat = NULL;
    allocate_matrix(&new_mat, self->mat->rows, self->mat->cols);
    result->mat = new_mat;
    result->shape = get_shape(self->mat->rows, self->mat->cols);

    int retval = abs_matrix(result->mat, self->mat);
    if (retval != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix absolute value failed");
        return NULL;
    }

    return (PyObject *) result;
}

/*
 * Raise numc.Matrix (Matrix61c) to the `pow`th power. You can ignore the argument `optional`.
 */
PyObject *Matrix61c_pow(Matrix61c *self, PyObject *pow, PyObject *optional) {

    if (!PyLong_Check(pow)) {
        PyErr_SetString(PyExc_TypeError, "Power argument must of type PyLong?");
        return NULL;
    }
    long p = PyLong_AsLong(pow);
    if (p < 0) {
        PyErr_SetString(PyExc_TypeError, "Power argument must be non-negative!");
        return NULL;
    }
    if (self->mat->rows != self->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Pow must be called on a square matrix!");
        return NULL;
    }

    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    matrix *new_mat = NULL;
    allocate_matrix(&new_mat, self->mat->rows, self->mat->cols);
    result->mat = new_mat;
    result->shape = get_shape(self->mat->rows, self->mat->cols);

    int retval = pow_matrix(result->mat, self->mat, p);
    if (retval != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Matrix power failed");
        return NULL;
    }

    return (PyObject *) result;
}

/*
 * Create a PyNumberMethods struct for overloading operators with all the number methods you have
 * defined. You might find this link helpful: https://docs.python.org/3.6/c-api/typeobj.html
 */
PyNumberMethods Matrix61c_as_number = {
    .nb_add = (binaryfunc) Matrix61c_add,
    .nb_subtract = (binaryfunc) Matrix61c_sub,
    .nb_multiply = (binaryfunc) Matrix61c_multiply,
    .nb_power = (ternaryfunc) Matrix61c_pow,
    .nb_absolute = (unaryfunc) Matrix61c_abs,
    .nb_negative = (unaryfunc) Matrix61c_neg
};


/* INSTANCE METHODS */

/*
 * Given a numc.Matrix self, parse `args` to (int) row, (int) col, and (double/int) val.
 * Return None in Python (this is different from returning null).
 */
PyObject *Matrix61c_set_value(Matrix61c *self, PyObject* args) {
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    if (PyArg_UnpackTuple(args, "args", 3, 3, &arg1, &arg2, &arg3)) {
        /* arguments are (rows, cols, val) */
        if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && (PyLong_Check(arg3)
                || PyFloat_Check(arg3))) {
	    if (PyLong_AsLong(arg1) < 0 || PyLong_AsLong(arg1) >= (self->mat->rows)) {
		PyErr_SetString(PyExc_IndexError, "Invalid row index");
	    } else if (PyLong_AsLong(arg2) < 0 || PyLong_AsLong(arg2) >= (self->mat->cols)) {
		PyErr_SetString(PyExc_IndexError, "Invalid col index");
	    } else {
		if (PyLong_Check(arg3)) {
                    (self->mat->data)[PyLong_AsLong(arg1)][PyLong_AsLong(arg2)] =  PyLong_AsLong(arg3);
		} else {
		    (self->mat->data)[PyLong_AsLong(arg1)][PyLong_AsLong(arg2)] =  PyFloat_AsDouble(arg3);
		}
	    }
	} else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
    }
    return Py_None;
}

/*
 * Given a numc.Matrix `self`, parse `args` to (int) row and (int) col.
 * Return the value at the `row`th row and `col`th column, which is a Python
 * float/int.
 */
PyObject *Matrix61c_get_value(Matrix61c *self, PyObject* args) {
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    if (PyArg_UnpackTuple(args, "args", 2, 2, &arg1, &arg2)) {
        /* arguments are (rows, cols) */
        if (arg1 && arg2 && PyLong_Check(arg1) && PyLong_Check(arg2)) {
            if (PyLong_AsLong(arg1) < 0 || PyLong_AsLong(arg1) >= (self->mat->rows)) {
                PyErr_SetString(PyExc_IndexError, "Invalid row index");
            } else if (PyLong_AsLong(arg2) < 0 || PyLong_AsLong(arg2) >= (self->mat->cols)) {
                PyErr_SetString(PyExc_IndexError, "Invalid col index");
            } else {
                return PyFloat_FromDouble((self->mat->data)[PyLong_AsLong(arg1)][PyLong_AsLong(arg2)]);
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
    }
    return Py_None;

}

/*
 * Create an array of PyMethodDef structs to hold the instance methods.
 * Name the python function corresponding to Matrix61c_get_value as "get" and Matrix61c_set_value
 * as "set"
 * You might find this link helpful: https://docs.python.org/3.6/c-api/structures.html
 */
PyMethodDef Matrix61c_methods[] = {
    /* TODO: YOUR CODE HERE */
    {"get", (PyCFunction) Matrix61c_get_value, METH_VARARGS, "get the value in the ith row and jth column"},
    {"set", (PyCFunction) Matrix61c_set_value, METH_VARARGS, "set the value in the ith row and jth column"},
    {NULL, NULL, 0, NULL}
};

/* INDEXING */

/*
 * Given a numc.Matrix `self`, index into it with `key`. Return the indexed result.
 */
PyObject *Matrix61c_subscript(Matrix61c* self, PyObject* key) {
    PyObject *key1 = NULL;
    PyObject *key2 = NULL;
    if (PyArg_UnpackTuple(key, "key", 2, 2, &key1, &key2)) {
	// In this case, we handle a tuple of ints and slices
    	if (self->mat->is_1d) {
	    // No 1D matrix should be indexed by a tuple
	    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
	    return Py_None;
	}
	if (PyLong_Check(key1) && PyLong_Check(key2)) {
	    // In this case, we handle [int, int]
	    return PyFloat_FromDouble((self->mat->data)[PyLong_AsLong(key1)][PyLong_AsLong(key2)]);

	} else if (PyLong_Check(key1) && PySlice_Check(key2)) {
	    // In this case, we handle [int, slice]
	    Py_ssize_t len = self->mat->cols, start1 = 0, stop1 = 0, step1 = 0, slicelength1 = 0;
	    PySlice_GetIndicesEx((PySliceObject*)key2,len,&start1,&stop1,&step1,&slicelength1);
	    long k = PyLong_AsLong(key1);
	    if (!valid_stepsize_and_slicelength(step1, slicelength1)) {
		// step1 must be 1 and slicelength1 must be greater than 0 
		PyErr_SetString(PyExc_ValueError, "Invalid arguments");
	    	return Py_None;
	    }
	    if (!valid_indices(self->mat->cols, start1, stop1)) {
		// the start and stop indices of the slice must be in bounds
		PyErr_SetString(PyExc_IndexError, "Index out of bounds");
		return Py_None;
	    }
	    if (!(k >= 0) || !(k < self->mat->rows)) {
		// the index must be in bounds
		PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return Py_None;
	    }
	    if (slicelength1 == 1) {
		// return a number, not a matrix
		return PyFloat_FromDouble((self->mat->data)[k][start1]);
	    }
	    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            matrix *new_mat = NULL;
            allocate_matrix_ref(&new_mat, self->mat, k, start1, 1, slicelength1);
            result->mat = new_mat;
            result->shape = get_shape(1, slicelength1);
            PyErr_SetString(NULL, "");
            return result;
	} else if (PySlice_Check(key1) && PyLong_Check(key2)) {
	    // In this case, we handle [slice, int]
	    Py_ssize_t len = self->mat->rows, start1 = 0, stop1 = 0, step1 = 0, slicelength1 = 0;
            PySlice_GetIndicesEx((PySliceObject*)key1,len,&start1,&stop1,&step1,&slicelength1);
	    long k = PyLong_AsLong(key2);
	    if (!valid_stepsize_and_slicelength(step1, slicelength1)) {
                // step1 must be 1 and slicelength1 must be greater than 0
	        PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return Py_None;
	    }
	    if (!valid_indices(self->mat->rows, start1, stop1)) {
		// the start and stop indices of the slice must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return Py_None;
	    }
	    if (!(k >= 0) || !(k < self->mat->cols)) {
                // the index must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return Py_None;
            }
	    if (slicelength1 == 1) {
                // return a number, not a matrix
                return PyFloat_FromDouble((self->mat->data)[start1][k]);
            }   
	    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            matrix *new_mat = NULL;
            allocate_matrix_ref(&new_mat, self->mat, start1, k, slicelength1, 1);
            result->mat = new_mat;
            result->shape = get_shape(slicelength1, 1);
            PyErr_SetString(NULL, "");
	    return result;
	} else if (PySlice_Check(key1) && PySlice_Check(key2)) {
	    // In this case, we handle [slice, slice]
	    Py_ssize_t len_r = self->mat->rows, start1_r = 0, stop1_r = 0, step1_r = 0, slicelength1_r = 0;
	    Py_ssize_t len_c = self->mat->cols, start1_c = 0, stop1_c = 0, step1_c = 0, slicelength1_c = 0;
	    PySlice_GetIndicesEx((PySliceObject*)key1,len_r,&start1_r,&stop1_r,&step1_r,&slicelength1_r);
	    PySlice_GetIndicesEx((PySliceObject*)key2,len_c,&start1_c,&stop1_c,&step1_c,&slicelength1_c);
	    if (!valid_stepsize_and_slicelength(step1_r, slicelength1_r) 
			    || !valid_stepsize_and_slicelength(step1_c, slicelength1_c)) {
                // step1 must be 1 and slicelength1 must be greater than 0
		PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return Py_None;
	    }
	    if (!valid_indices(self->mat->cols, start1_c, stop1_c) || !valid_indices(self->mat->rows, start1_r, stop1_r)) {
                // the start and stop indices of the slices must be in bounds
		PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return Py_None;
	    }
	    if (slicelength1_r == 1 && slicelength1_c == 1) {
                // return a number, not a matrix
                return PyFloat_FromDouble((self->mat->data)[start1_r][start1_c]);
            }
            Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            matrix *new_mat = NULL;
            allocate_matrix_ref(&new_mat, self->mat, start1_r, start1_c, slicelength1_r, slicelength1_c);
            result->mat = new_mat;
            result->shape = get_shape(slicelength1_r, slicelength1_c);
            PyErr_SetString(NULL, "");
	    return result;
	} else {
	    // if the tuple does not parse to ints and slices, throw an error
	    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
	    return Py_None;
	}
    } else if (PyLong_Check(key)) {
	// In this case, we handle [int]
	PyErr_SetString(NULL, ""); // this resets the PyErr to be empty (needed)
	long k = PyLong_AsLong(key); 
	if (self->mat->is_1d) {
	    // Here we handle the 1D case
            if (!(k >= 0) || !(k < self->mat->rows * self->mat->cols)) {
		// the index must be in bounds
		PyErr_SetString(PyExc_IndexError, "Index out of bounds");
	   	return Py_None;
	    }
	    if (self->mat->rows == 1) {
		// Here we handle a matrix that looks like a row vector
                return PyFloat_FromDouble((self->mat->data)[0][k]);
	    } else {
		// Here we handle a matrix that looks like a column vector
		return PyFloat_FromDouble((self->mat->data)[k][0]);
	    }
	} else {
            if (!(k >= 0) || !(k < self->mat->rows)) {
                // the index must be in bounds
		PyErr_SetString(PyExc_IndexError, "Index out of bounds");
            	return Py_None;
	    }
            Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            matrix *new_mat = NULL;
            allocate_matrix_ref(&new_mat, self->mat, k, 0, 1, self->mat->cols);
            result->mat = new_mat;
            result->shape = get_shape(1, self->mat->cols);
            PyErr_SetString(NULL, "");
            return result;
	}

    } else if (PySlice_Check(key)) {
	// In this case, we handle [slice]
	PyErr_SetString(NULL, ""); // this resets the PyErr to be empty (needed)
	Py_ssize_t len = self->mat->rows, start1 = 0, stop1 = 0, step1 = 0, slicelength1 = 0;
        
	if (self->mat->is_1d) {
	    len = self->mat->rows * self->mat->cols;
	}
	
	PySlice_GetIndicesEx((PySliceObject*)key,len,&start1,&stop1,&step1,&slicelength1);
	if (!valid_stepsize_and_slicelength(step1, slicelength1)) {
		// step1 must be 1 and slicelength1 must be greater than 0
		PyErr_SetString(PyExc_ValueError, "Invalid arguments");
        	return Py_None;
	}

	if (self->mat->is_1d) {
	    // Here, we handle the 1D case
	    if (!valid_indices(self->mat->rows * self->mat->cols, start1, stop1)) {
		// the start and stop indices of the slices must be in bounds
		PyErr_SetString(PyExc_IndexError, "Index out of bounds");
	    	return Py_None;
	    }
	    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            matrix *new_mat = NULL;
	    if (self->mat->rows == 1) {
		// Here we handle a matrix that looks like a row vector
		if (slicelength1 == 1) {
		    return PyFloat_FromDouble((self->mat->data)[0][start1]);
		}
	        allocate_matrix_ref(&new_mat, self->mat, 0, start1, self->mat->rows, slicelength1);
	    } else {
		// Here we handle a matrix that looks like a column vector
		if (slicelength1 == 1) {
                    return PyFloat_FromDouble((self->mat->data)[start1][0]);
                }
		allocate_matrix_ref(&new_mat, self->mat, start1, 0, slicelength1, self->mat->cols);
	    }
	    result->mat = new_mat;
            result->shape = get_shape(slicelength1, self->mat->cols);
	    
	    return result;
	} else {
	    // Here, we handle the 2D case
	    if (!valid_indices(self->mat->rows, start1, stop1)) {
                // the start and stop indices of the slices must be in bounds
		PyErr_SetString(PyExc_IndexError, "Index out of bounds");
		return Py_None;
            }
	    Matrix61c *result = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
	    matrix *new_mat = NULL;
   	    allocate_matrix_ref(&new_mat, self->mat, start1, 0, slicelength1, self->mat->cols);
	    result->mat = new_mat;
   	    result->shape = get_shape(slicelength1, self->mat->cols);
	    PyErr_SetString(NULL, "");
	    return result;
	}
    } else {
	// Reaching this case means that args did not match any of the expected formats
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return Py_None;
    }
    return Py_None;
}

int valid_stepsize_and_slicelength(int step, int slicelength) {
    int valid = 1;
    if (step != 1) {
	valid = 0;
    }
    if (slicelength < 1) {
	valid = 0;
    }
    return valid;
}   

int valid_indices(int upper, int start, int stop) {
    int valid = 1;
    if (!(0 <= start && start <= upper) || !(0 <= stop && stop <= upper)) {
        valid = 0;
    }
    return valid;
}

/*
 * Given a numc.Matrix `self`, index into it with `key`, and set the indexed result to `v`.
 */
int Matrix61c_set_subscript(Matrix61c* self, PyObject *key, PyObject *v) {
    PyObject *key1 = NULL;
    PyObject *key2 = NULL;
    if (PyArg_UnpackTuple(key, "key", 2, 2, &key1, &key2)) {
        // In this case, we handle a tuple of ints and slices
        if (self->mat->is_1d) {
            // No 1D matrix should be indexed by a tuple
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (PyLong_Check(key1) && PyLong_Check(key2)) {
            // In this case, we handle [int, int]
	    if (PyLong_Check(v) || PyFloat_Check(v)) {
		(self->mat->data)[PyLong_AsLong(key1)][PyLong_AsLong(key2)] = PyFloat_AsDouble(v);
	    	return 0;
	    } else {
		PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
	    }
	} else if (PyLong_Check(key1) && PySlice_Check(key2)) {
            // In this case, we handle [int, slice]
            Py_ssize_t len = self->mat->cols, start1 = 0, stop1 = 0, step1 = 0, slicelength1 = 0;
            PySlice_GetIndicesEx((PySliceObject*)key2,len,&start1,&stop1,&step1,&slicelength1);
            long k = PyLong_AsLong(key1);
            if (!valid_stepsize_and_slicelength(step1, slicelength1)) {
                // step1 must be 1 and slicelength1 must be greater than 0 
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
            }
            if (!valid_indices(self->mat->cols, start1, stop1)) {
                // the start and stop indices of the slice must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return -1;
            }
            if (!(k >= 0) || !(k < self->mat->rows)) {
                // the index must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return -1;
            }
            if (slicelength1 == 1) {
                // return a number, not a matrix
		if (PyLong_Check(v) || PyFloat_Check(v)) {
                    (self->mat->data)[k][start1] =  PyFloat_AsDouble(v);
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                    return -1;
                }
            }
	    if (!PyList_Check(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
            if (PyList_Size(v) != slicelength1) {
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
            }
            if (!all_elements_int_float(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
            for (int i = start1; i < stop1; i++) {
                (self->mat->data)[k][i] = PyFloat_AsDouble(PyList_GetItem(v, i-start1));
            }
            return 0;

        } else if (PySlice_Check(key1) && PyLong_Check(key2)) {
            // In this case, we handle [slice, int]
            Py_ssize_t len = self->mat->rows, start1 = 0, stop1 = 0, step1 = 0, slicelength1 = 0;
            PySlice_GetIndicesEx((PySliceObject*)key1,len,&start1,&stop1,&step1,&slicelength1);
            long k = PyLong_AsLong(key2);
            if (!valid_stepsize_and_slicelength(step1, slicelength1)) {
                // step1 must be 1 and slicelength1 must be greater than 0
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
		return -1;
            }
            if (!valid_indices(self->mat->rows, start1, stop1)) {
                // the start and stop indices of the slice must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return -1;
	    }
            if (!(k >= 0) || !(k < self->mat->cols)) {
                // the index must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return -1;
            }
            if (slicelength1 == 1) {
                // return a number, not a matrix
		if (PyLong_Check(v) || PyFloat_Check(v)) {
                    (self->mat->data)[start1][k] =  PyFloat_AsDouble(v);
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                    return -1;
                }
            }
	    if (!PyList_Check(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
            if (PyList_Size(v) != slicelength1) {
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
            }
            if (!all_elements_int_float(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
            for (int i = start1; i < stop1; i++) {
                (self->mat->data)[i][k] = PyFloat_AsDouble(PyList_GetItem(v, i-start1));
            }
            return 0;

        } else if (PySlice_Check(key1) && PySlice_Check(key2)) {
            // In this case, we handle [slice, slice]
            Py_ssize_t len_r = self->mat->rows, start1_r = 0, stop1_r = 0, step1_r = 0, slicelength1_r = 0;
            Py_ssize_t len_c = self->mat->cols, start1_c = 0, stop1_c = 0, step1_c = 0, slicelength1_c = 0;
            PySlice_GetIndicesEx((PySliceObject*)key1,len_r,&start1_r,&stop1_r,&step1_r,&slicelength1_r);
            PySlice_GetIndicesEx((PySliceObject*)key2,len_c,&start1_c,&stop1_c,&step1_c,&slicelength1_c);
            if (!valid_stepsize_and_slicelength(step1_r, slicelength1_r)
                            || !valid_stepsize_and_slicelength(step1_c, slicelength1_c)) {
                // step1 must be 1 and slicelength1 must be greater than 0
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
	    }
            if (!valid_indices(self->mat->cols, start1_c, stop1_c) || !valid_indices(self->mat->rows, start1_r, stop1_r)) {
                // the start and stop indices of the slices must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return -1;
	    }
            if (slicelength1_r == 1 && slicelength1_c == 1) {
                // return a number, not a matrix
		if (PyLong_Check(v) || PyFloat_Check(v)) {
                    (self->mat->data)[start1_r][start1_c] = PyFloat_AsDouble(v);
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                    return -1;
                }
            }
	    
	    if (!PyList_Check(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
	    if (slicelength1_r == 1) {
	        // handle one row, multiple columns 	
	    	if (PyList_Size(v) != slicelength1_c ||!all_elements_int_float(v)) {
                    PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                    return -1;
                }
		for (int j = start1_c; j < stop1_c; j++) {
                    (self->mat->data)[start1_r][j] = PyFloat_AsDouble(PyList_GetItem(v, j-start1_c));
                }
		return 0;
	    }

	    if (slicelength1_c == 1) {
	        // handle multiple rows, one column
		if (PyList_Size(v) != slicelength1_r ||!all_elements_int_float(v)) {
                    PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                    return -1;
                }
                for (int j = start1_r; j < stop1_r; j++) {
                    (self->mat->data)[j][start1_c] = PyFloat_AsDouble(PyList_GetItem(v, j-start1_r));
                }
                return 0;
	    }

            if (PyList_Size(v) != slicelength1_r) {
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
            }
            if (!all_elements_list(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
            for (int i=0; i < slicelength1_r; i++) {
                if (!(slicelength1_c  == PyList_Size((PyList_GetItem(v, i))))
                                        || !all_elements_int_float(PyList_GetItem(v, i))) {
                    PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                    return -1;
                }
            }
            for (int i = start1_r; i < stop1_r; i++) {
                for (int j = start1_c; j < stop1_c; j++) {
                    (self->mat->data)[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(v, i-start1_r), j-start1_c));
                }
            }
            return 0;
        } else {
            // if the tuple does not parse to ints and slices, throw an error
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
	    return -1;
        }
    } else if (PyLong_Check(key)) {
        // In this case, we handle [int]
        PyErr_SetString(NULL, ""); // this resets the PyErr to be empty (needed)
        long k = PyLong_AsLong(key);
	if (self->mat->is_1d) {
            // Here we handle the 1D case
            if (!(k >= 0) || !(k < self->mat->rows * self->mat->cols)) {
                // the index must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return -1;
            }
            if (self->mat->rows == 1) {
                // Here we handle a matrix that looks like a row vector
		if (PyLong_Check(v) || PyFloat_Check(v)) {
                    (self->mat->data)[0][k] = PyFloat_AsDouble(v);
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                    return -1;
                }

            } else {
                // Here we handle a matrix that looks like a column vector
		if (PyLong_Check(v) || PyFloat_Check(v)) {
                    (self->mat->data)[k][0] = PyFloat_AsDouble(v);
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                    return -1;
                }
            }
        } else {
	    // Here we handle the 2D case
            if (!(k >= 0) || !(k < self->mat->rows)) {
                // the index must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return -1;
            }
	    if (!PyList_Check(v)) {
		PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
	    }
	    if (PyList_Size(v) != self->mat->cols) {
		PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
	    }
	    if (!all_elements_int_float(v)) {
		PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
	    }
	    int x = 0; // x=1 for longs, x=0 for floats
	    if (PyLong_Check(PyList_GetItem(v, 0))) {
		x = 1;
	    }
	    for (int i = 0; i < self->mat->cols; i++) {
		(self->mat->data)[k][i] = PyFloat_AsDouble(PyList_GetItem(v, i));
	    }
	    return 0;
        }

    } else if (PySlice_Check(key)) {
        // In this case, we handle [slice]
        PyErr_SetString(NULL, ""); // this resets the PyErr to be empty (needed)
        Py_ssize_t len = self->mat->rows, start1 = 0, stop1 = 0, step1 = 0, slicelength1 = 0;

        if (self->mat->is_1d) {
            len = self->mat->rows * self->mat->cols;
        }

        PySlice_GetIndicesEx((PySliceObject*)key,len,&start1,&stop1,&step1,&slicelength1);
        if (!valid_stepsize_and_slicelength(step1, slicelength1)) {
                // step1 must be 1 and slicelength1 must be greater than 0
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
        }

        if (self->mat->is_1d) {
            // Here, we handle the 1D case
            if (!valid_indices(self->mat->rows * self->mat->cols, start1, stop1)) {
                // the start and stop indices of the slices must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds");
                return -1;
            }

	    if (slicelength1 == 1 && self->mat->rows==1) {
		// handle the slice of one number for column vector case
                if (PyLong_Check(v) || PyFloat_Check(v)) {
                    (self->mat->data)[0][start1] = PyFloat_AsDouble(v);
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                    return -1;
                }
            }
	    if (slicelength1 == 1 && self->mat->cols==1) {
                if (PyLong_Check(v) || PyFloat_Check(v)) {
                    (self->mat->data)[start1][0] = PyFloat_AsDouble(v);
                    return 0;
                } else {
                    PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                    return -1;
                }
            }
            
            if (!PyList_Check(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
            if (PyList_Size(v) != slicelength1) {
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
            }
            if (!all_elements_int_float(v)) {
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
            }

	    if (self->mat->rows == 1) {
                // Here we handle a matrix that looks like a row vector
		for (int i = start1; i < stop1; i++) {
                    (self->mat->data)[0][i] = PyFloat_AsDouble(PyList_GetItem(v, i-start1));
                }
		return 0;
	    } else {
                // Here we handle a matrix that looks like a column vector
		for (int i = start1; i < stop1; i++) {
                    (self->mat->data)[i][0] = PyFloat_AsDouble(PyList_GetItem(v, i-start1));
                }
                return 0;
	    }
        } else {
            // Here, we handle the 2D case
            if (!valid_indices(self->mat->rows, start1, stop1)) {
                // the start and stop indices of the slices must be in bounds
                PyErr_SetString(PyExc_IndexError, "Index out of bounds \n");
                return -1;
            }
	    if (!PyList_Check(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
            if (PyList_Size(v) != slicelength1) {
                PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                return -1;
            }
            if (!all_elements_list(v)) {
                PyErr_SetString(PyExc_TypeError, "Invalid arguments");
                return -1;
            }
	    for (int i=0; i < slicelength1; i++) {
		if (!(self->mat->cols == PyList_Size((PyList_GetItem(v, i)))) 
					|| !all_elements_int_float(PyList_GetItem(v, i))) {
		    PyErr_SetString(PyExc_ValueError, "Invalid arguments");
                    return -1;
		}
	    }
	    for (int i = start1; i < stop1; i++) {
		for (int j = 0; j < self->mat->cols; j++) {
		    (self->mat->data)[i][j] = PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(v, i-start1), j));	    
		}
	    }
	    return 0;
	}
    } else {
        // Reaching this case means that args did not match any of the expected formats
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return -1;
    }
    return -1;
}

int all_elements_int_float(PyObject *ls) {
    int x = 1;
    for (int i=0; i < PyList_Size(ls); i++) {
        if (!(PyLong_Check(PyList_GetItem(ls, i))) && !(PyFloat_Check(PyList_GetItem(ls, i)))) {
	    x = 0;
        }
    }
    return x;
}

int all_elements_list(PyObject *ls) {
    int x = 1;
    for (int i=0; i < PyList_Size(ls); i++) {
	if (!PyList_Check(PyList_GetItem(ls, i))) {
	    x = 0;
	}
    }
    return x;
}
    
PyMappingMethods Matrix61c_mapping = {
    NULL,
    (binaryfunc) Matrix61c_subscript,
    (objobjargproc) Matrix61c_set_subscript,
};

/* INSTANCE ATTRIBUTES*/
PyMemberDef Matrix61c_members[] = {
    {
        "shape", T_OBJECT_EX, offsetof(Matrix61c, shape), 0,
        "(rows, cols)"
    },
    {NULL}  /* Sentinel */
};

PyTypeObject Matrix61cType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "numc.Matrix",
    .tp_basicsize = sizeof(Matrix61c),
    .tp_dealloc = (destructor)Matrix61c_dealloc,
    .tp_repr = (reprfunc)Matrix61c_repr,
    .tp_as_number = &Matrix61c_as_number,
    .tp_flags = Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,
    .tp_doc = "numc.Matrix objects",
    .tp_methods = Matrix61c_methods,
    .tp_members = Matrix61c_members,
    .tp_as_mapping = &Matrix61c_mapping,
    .tp_init = (initproc)Matrix61c_init,
    .tp_new = Matrix61c_new
};


struct PyModuleDef numcmodule = {
    PyModuleDef_HEAD_INIT,
    "numc",
    "Numc matrix operations",
    -1,
    Matrix61c_class_methods
};

/* Initialize the numc module */
PyMODINIT_FUNC PyInit_numc(void) {
    PyObject* m;

    if (PyType_Ready(&Matrix61cType) < 0)
        return NULL;

    m = PyModule_Create(&numcmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&Matrix61cType);
    PyModule_AddObject(m, "Matrix", (PyObject *)&Matrix61cType);
    printf("CS61C Fall 2020 Project 4: numc imported!\n");
    fflush(stdout);
    return m;
}
