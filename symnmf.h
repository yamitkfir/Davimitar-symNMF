#ifndef SYMNMF_H
#define SYMNMF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symnmf.c"

/* Function declarations */
double squared_euclidean_dist(double* point1, double* point2, int dimension);
double** get_column(double** M, int rows_num, int j);
double** optimizing_H(double** H, int rows_num, int cols_num, double** W);
int update_H(double** W, double** H, double** new_H, int n, int k);
double** similarity_matrix(double** datapoints, int n, int d);
double** diagonal_degree_matrix(double** A, int n);
double** normalized_similarity_matrix(double** sim_matrix, int n);

double **read_data(const char *filename, int *n, int *d);
void print_matrix(double **matrix, int rows, int cols);

/* Helper functions */
void run_selected_algorithm(const char* goal, double** A, double** points, int n, int d, double** result_pointer);
double **create_points_matrix(FILE *fp, char line[], int *n, int *d);
double sq_frobenius_norm(double** A, int rows_num, int cols_num, double** B);
void free_matrix(double** M, int len);
double** multiply_matrix(double** matrixA, double** matrixB, int m, int n, int k);
double matrix_mult_cell(double** A, int A_cols_num ,double** B, int i, int j);
void update_H_cell(double **W, double **H, double **new_H, double **HtH_col, int n, int k, int i, int j);
void exit_with_error();

#endif