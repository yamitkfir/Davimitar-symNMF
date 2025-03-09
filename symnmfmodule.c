/* Python C API Wrapper
In this file you will define your C extension which will serve the functions:
symnmf,sym,ddg,norm for Python */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmfalgo.h"

/* Function declarations */



static PyObject* symnmf(PyObject* self, PyObject* args) {
    // TODO David
    printf("WIP\n");
    Py_RETURN_NONE;
}

static PyObject* sym(PyObject* self, PyObject* args) {
    // TODO David
    printf("WIP\n");
    Py_RETURN_NONE;
}

static PyObject* ddg(PyObject* self, PyObject* args) {
    // TODO David
    printf("WIP\n");
    Py_RETURN_NONE;
}

static PyObject* norm(PyObject* self, PyObject* args) {
    // TODO David
    printf("WIP\n");
    Py_RETURN_NONE;
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
