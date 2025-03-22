/*
Python C API Wrapper
In this file you will define your C extension which will serve the functions:
symnmf,sym,ddg,norm for Python
*/

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

#define ERR_LIST_FORMAT "Expected a list of lists of floats"
#define ERR_LIST_ITEM_FORMAT "List items must be floats"
#define ERR_SYMNMF_FORMAT "Input must be two matrixes"

/* Function declarations - for module use only */
static PyObject* symnmf(PyObject* self, PyObject* args);
static PyObject* sym(PyObject* self, PyObject* args);
static PyObject* ddg(PyObject* self, PyObject* args);
static PyObject* norm(PyObject* self, PyObject* args);
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
    PyObject *lstH, *lstW, *ret;
    double** H, **W;
    int n, k;
    printf("Enter"); /* TODO saar pls dont forget to delete these prints when youre done with them */
    if(!PyArg_ParseTuple(args, "OO", &lstW, &lstH)) {
        PyErr_SetString(PyExc_TypeError, ERR_SYMNMF_FORMAT);
        Py_RETURN_NONE;
    }
    if (!PyList_Check(lstH) || !PyList_Check(lstW)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    printf("Legal");
    H = getDataPoints(lstH);
    printf("H");
    W = getDataPoints(lstW);
    printf("W");
    if(H == NULL || W == NULL) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    n = PyList_Size(lstH);
    k = PyList_Size(PyList_GetItem(lstH, 0));
    printf("enter opti");
    optimizing_H(H, n, k, W); 
    printf("Got H");
    freeDataPoints(W, n);
    printf("free W");
    ret = MatrixToPyList(H, n, k);
    freeDataPoints(H, n);
    printf("free H");
    return ret;
}

/*
Input: Datapoints Py List
Output: Similarity matrix
Given a matrix of datapoints, calculate the similarity matrix.
Stage 1.1 in the instructions.
*/
static PyObject* sym(PyObject* self, PyObject* args) {
    PyObject* lst, *ret;
    double** A, **dataPoints;
    if(!PyArg_ParseTuple(args, "O", &lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    if (!PyList_Check(lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    dataPoints = getDataPoints(lst);
    if(dataPoints == NULL) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    A = similarity_matrix(dataPoints, PyList_Size(lst), PyList_Size(PyList_GetItem(lst, 0)));
    freeDataPoints(dataPoints, PyList_Size(lst));
    
    ret = MatrixToPyList(A, PyList_Size(lst), PyList_Size(lst));
    free_matrix(A, PyList_Size(lst));
    return ret;
}

/*
Input: Datapoints Py List
Output: Diagonal Degree Matrix
Given a Datapoints matrix, calculate the Diagonal Degree Matrix.
Stage 1.2 in the instructions.
*/
static PyObject* ddg(PyObject* self, PyObject* args) {
    PyObject* lst, *ret;
    double** dataPoints, **D, **A;
    if(!PyArg_ParseTuple(args, "O", &lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    if (!PyList_Check(lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    dataPoints = getDataPoints(lst);
    if(dataPoints == NULL) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    A = similarity_matrix(dataPoints, PyList_Size(lst), PyList_Size(PyList_GetItem(lst, 0)));
    D = diagonal_degree_matrix(A, PyList_Size(lst));
    freeDataPoints(dataPoints, PyList_Size(lst));
    free_matrix(A, PyList_Size(lst));
    ret = MatrixToPyList(D, PyList_Size(lst), PyList_Size(lst));
    free_matrix(D, PyList_Size(lst));
    return ret;
}

/*
Input: Datapoints Py List
Output: Normalized Similarity Matrix
Given a Datapoints matrix, calculate the Normalized Similarity Matrix.
Stage 1.3 in the instructions.
*/
static PyObject* norm(PyObject* self, PyObject* args) {
    PyObject* lst, *ret;
    double** dataPoints, **A, **normalized;
    if(!PyArg_ParseTuple(args, "O", &lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    if (!PyList_Check(lst)) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    dataPoints = getDataPoints(lst);
    if(dataPoints == NULL) {
        PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
        Py_RETURN_NONE;
    }
    A = similarity_matrix(dataPoints, PyList_Size(lst), PyList_Size(PyList_GetItem(lst, 0)));
    normalized = normalized_similarity_matrix(A, PyList_Size(lst));
    freeDataPoints(dataPoints, PyList_Size(lst));
    free_matrix(A, PyList_Size(lst));
    ret = MatrixToPyList(normalized, PyList_Size(lst), PyList_Size(lst));

    free_matrix(normalized, PyList_Size(lst));
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
    int i;
    for (i = 0; i < n; i++) {
        free(dataPoints[i]);
    }
    free(dataPoints);
}

double** getDataPoints(PyObject* lst) {
    Py_ssize_t len = PyList_Size(lst);
    double** dataPoints = malloc(len * sizeof(double*)); /* TODO david replaced calloc with malloc, ok? */
    Py_ssize_t i, j;
    for (i = 0; i < len; i++) {
        PyObject* subList = PyList_GetItem(lst, i);
        if (!PyList_Check(subList)) {
            PyErr_SetString(PyExc_TypeError, ERR_LIST_FORMAT);
            return NULL;
        }
        Py_ssize_t subListLen = PyList_Size(subList);
        dataPoints[i] = malloc(subListLen * sizeof(double)); /* TODO david replaced calloc with malloc, ok? */
        for(j = 0; j < subListLen; j++) {
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
    PyObject* lst = PyList_New(n), *num, *subList;
    int i, j;
    for (i = 0; i < n; i++) {
        subList = PyList_New(m);
        for (j = 0; j < m; j++) {
            num = PyFloat_FromDouble(matrix[i][j]);
            PyList_SET_ITEM(subList, j, num);
        }
        PyList_SET_ITEM(lst, i, subList);
    }
    return lst;
}
