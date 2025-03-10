/* Python C API Wrapper
In this file you will define your C extension which will serve the functions:
symnmf,sym,ddg,norm for Python */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmfalgo.h"

#define ERR_LIST_FORMAT "Expected a list of lists of floats"
#define ERR_LIST_ITEM_FORMAT "List items must be floats"
#define ERR_SYMNF_FORMAT "Input must be two matrixes"

/* Function declarations - for module use only */
double** getDataPoints(PyObject* lst);
void freeDataPoints(double** dataPoints, int n);
PyObject* MatrixToPyList(double** matrix, int n, int m);

/*
Input: Matrixes H, W and more?
Output: Final H
Given a starting matrix H and a graph laplacian W, perform the optimization algorithm in the instructions.
Stages 1.4 and 1.5 in the instructions.
*/
static PyObject* symnmf(PyObject* self, PyObject* args) {
    // TODO David
    printf("WIP\n");
    PyObject* lst;
    if(!PyArg_ParseTuple(args, "O", &lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_SYMNF_FORMAT);
        Py_RETURN_NONE;
    }
    if (!PyList_Check(lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_SYMNF_FORMAT);
        Py_RETURN_NONE;
    }
    if(PyList_Size(lst) != 2) {
        PyErr_SetString(PyExc_TypeError, ERR_SYMNF_FORMAT);
        Py_RETURN_NONE;
    }
    if(!PyList_Check(PyList_GetItem(lst, 0)) || !PyList_Check(PyList_GetItem(lst, 1))) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    double** H = getDataPoints(PyList_GetItem(lst, 0));
    double** W = getDataPoints(PyList_GetItem(lst, 1));
    if(H == NULL || W == NULL) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    int n = PyList_Size(PyList_GetItem(lst, 0));
    int k = PyList_Size(PyList_GetItem(PyList_GetItem(lst, 0, 0)));
    double** new_H = optimizing_H(H, n, k, W); // waiting for Response to determine rows v cols
    freeDataPoints(H, n);
    freeDataPoints(W, n);
    PyObject* ret = MatrixToPyList(new_H, n, k);
    free_matrix(new_H);
    return ret;
}

/*
Input: Datapoints Py List
Output: Similarity matrix
Given a matrix of datapoints, calculate the similarity matrix.
Stage 1.1 in the instructions.
*/
static PyObject* sym(PyObject* self, PyObject* args) {
    printf("WIP\n"); /* To Delete */
    PyObject* lst;
    if(!PyArg_ParseTuple(args, "O", &lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    if (!PyList_Check(lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    double** dataPoints = getDataPoints(lst);
    if(dataPoints == NULL) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    double** A = similarity_matrix(dataPoints, PyList_Size(lst), PyList_Size(PyList_GetItem(lst, 0)));
    freeDataPoints(dataPoints, PyList_Size(lst));
    
    PyObject* ret = MatrixToPyList(A, PyList_Size(lst), PyList_Size(lst));
    free_matrix(A);
    return ret;
}

/*
Input: Datapoints Py List
Output: Diagonal Degree Matrix
Given a Datapoints matrix, calculate the Diagonal Degree Matrix.
Stage 1.2 in the instructions.
*/
static PyObject* ddg(PyObject* self, PyObject* args) {
    printf("WIP\n");
    PyObject* lst;
    if(!PyArg_ParseTuple(args, "O", &lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    if (!PyList_Check(lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    double** dataPoints = getDataPoints(lst);
    if(dataPoints == NULL) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    double** A = similarity_matrix(dataPoints, PyList_Size(lst), PyList_Size(PyList_GetItem(lst, 0)));
    double** D = diagonal_degree_matrix(A, PyList_Size(lst));
    freeDataPoints(dataPoints, PyList_Size(lst));
    free_matrix(A);
    PyObject* ret = MatrixToPyList(D, PyList_Size(lst), PyList_Size(lst));
    free_matrix(D);
    return ret;
}

/*
Input: Datapoints Py List
Output: Normalized Similarity Matrix
Given a Datapoints matrix, calculate the Normalized Similarity Matrix.
Stage 1.3 in the instructions.
*/
static PyObject* norm(PyObject* self, PyObject* args) {
    printf("WIP\n");
    PyObject* lst;
    if(!PyArg_ParseTuple(args, "O", &lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    if (!PyList_Check(lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    double** dataPoints = getDataPoints(lst);
    if(dataPoints == NULL) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    double** A = similarity_matrix(dataPoints, PyList_Size(lst), PyList_Size(PyList_GetItem(lst, 0)));
    double** normalized = normalized_similarity_matrix(A, PyList_Size(lst));
    freeDataPoints(dataPoints, PyList_Size(lst));
    free_matrix(A);
    PyObject* ret = MatrixToPyList(normalized, PyList_Size(lst), PyList_Size(lst));
    free_matrix(normalized);
    return ret;
}

static PyMethodDef symnmfmethods[] = {
    {"symnmf", symnmf, METH_VARARGS, "Performs SymNMF on a matrix."},
    {"sym", sym, METH_VARARGS, "Performs Sym on a matrix."},
    {"ddg", ddg, METH_VARARGS, "Performs DDG on a matrix."},
    {"norm", norm, METH_VARARGS, "Performs Norm on a matrix."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    "A module that performs SymNMF on a matrix.",
    -1,
    symnmfmethods
};

PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject* m = PyModule_Create(&symnmfmodule);
    if (m == NULL) {
        return NULL;
    }
    return m;
}

void freeDataPoints(double** dataPoints, int n) {
    for (int i = 0; i < n; i++) {
        free(dataPoints[i]);
    }
    free(dataPoints);
}

double** getDataPoints(PyObject* lst) {
    Py_ssize_t len = PyList_Size(lst);
    double** dataPoints = calloc(len, sizeof(double*));

    for (Py_ssize_t i = 0; i < len; i++) {
        PyObject* subList = PyList_GetItem(lst, i);
        if (!PyList_Check(subList)) {
            PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
            return NULL;
        }
        Py_ssize_t subListLen = PyList_Size(subList);
        dataPoints[i] = calloc(subListLen, sizeof(double));
        for(Py_ssize_t j = 0; j < subListLen; j++) {
            PyObject* cord = PyList_GetItem(subList, j);
            if (!PyFloat_Check(cord) && !PyLong_Check(cord)) {
                PyErr_SetString(PyExc_TypeError, ERR_LIST_ITEM_FORMAT);
                return NULL;
            }
            dataPoints[i][j] = PyFloat_AsDouble(cord);
        }
    }
    return dataPoints;
}

PyObject* MatrixToPyList(double** matrix, int n, int m) {
    PyObject* lst = PyList_New(n);
    for (int i = 0; i < n; i++) {
        PyObject* subList = PyList_New(m);
        for (int j = 0; j < m; j++) {
            PyList_SetItem(subList, j, PyFloat_FromDouble(matrix[i][j]));
        }
        PyList_SetItem(lst, i, subList);
    }
    return lst;
}
