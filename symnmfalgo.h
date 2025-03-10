#ifndef SYMNMFALGO_H
#define SYMNMFALGO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Function declarations */
double squared_euclidean_dist(double* point1, double* point2, int dimension);
double** get_column(double** M, int rows_num, int j);
double** optimizing_H(double** H, int rows_num, int cols_num, double** W);
int update_H(double** W, double** H, double** new_H, int n, int k);
double** similarity_matrix(double** datapoints, int n, int d);
double** diagonal_degree_matrix(double** A, int n);

/* Helper functions */
double sq_euclidian_norm(double* point1, double* point2, int dimension);
double frobenius_norm(double** A, int rows_num, int cols_num, double** B);
void free_matrix(double** M, int len);
double** multiplyMatrix(double** matrixA, double** matrixB, int m, int n, int k); /* A - m x n, B - n x k */
double matrix_mult_cell(double** A, int A_cols_num ,double** B, int i, int j);
void exit_with_error();

#endif /* SYMNMFALGO_H */