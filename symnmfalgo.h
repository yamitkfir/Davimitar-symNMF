#ifndef SYMNMFALGO_H
#define SYMNMFALGO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Function declarations */
void free_matrix(double** M);
void exit_with_error();
double matrix_mult_cell(double** A, double** B, int i, int j);
double squared_euclidean_dist(double* point1, double* point2, int dimension);
double frobenius_norm(double** A, double** B);
double** get_column(double** M, int j);
double** optimizing_H(double** H, double** W);
int update_H(double** W, double** H, double** new_H, int n, int k);
double** similarity_matrix(double** datapoints, int n, int d);
double** diagonal_degree_matrix(double** A, int n);
double** multiplyMatrix(double** matrixA, double** matrixB, int m, int n, int k); // A - m x n, B - n x k
double** similarity_matrix(double** datapoints, int n, int d);
double** normalized_similarity_matrix(double** sim_matrix, int n);

#endif /* SYMNMFALGO_H */