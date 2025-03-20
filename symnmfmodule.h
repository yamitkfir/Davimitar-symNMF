#ifndef SYMNMFMODULE_H
#define SYMNMFMODULE_H

#include <Python.h>
#include "symnmf.h"

/* Function declarations */
static PyObject* symnmf(PyObject* self, PyObject* args);
static PyObject* sym(PyObject* self, PyObject* args);
static PyObject* ddg(PyObject* self, PyObject* args);
static PyObject* norm(PyObject* self, PyObject* args);
double** getDataPoints(PyObject* lst);
void freeDataPoints(double** dataPoints, int n);
PyObject* MatrixToPyList(double** matrix, int n, int m);

#endif