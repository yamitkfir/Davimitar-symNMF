#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symnmfalgo.h"

#define max_iter 300
#define eps 1e-4
#define beta 0.5
#define ERROR_MSG "An Error Has Occurred\n"

/* 
TODO diagonal_matrix()
TODO David: normalized_similarity_matrix()
TODO init_H()
*/

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

/*
* Calculate the Diagonal Degree Matrix D for a given similarity matrix A.
* Uses a 2D array representation (array of arrays).
*/
double** diagonal_degree_matrix(double** A, int n) {
    double** D = (double**)malloc(n * sizeof(double*));
    int i, j; double sum;
    if (D == NULL) {
        exit_with_error();
    }
    /* Allocate memory for each row and init to 0 */
    for (i = 0; i < n; i++) {
        D[i] = (double*)calloc(n, sizeof(double));
        if (D[i] == NULL) {
            free_matrix(D, i);
            exit_with_error();
        }
    }
    /* Calc degrees */
    for (i = 0; i < n; i++) {
        sum = 0.0;
        /* Sum the i-th row of A to get the degree */
        for (j = 0; j < n; j++) {
            sum += A[i][j];
        }
        D[i][i] = sum; /* All other elements remain zero (from calloc) */
    }
    return D;
}

/* TODO free_matrix(A) once done with it! */
double** similarity_matrix(double** datapoints, int n, int d) {
    double** A = (double**)malloc(n * sizeof(double*));
    int i, j; double dist;
    if (A == NULL) {
        exit_with_error();
    }
    /* Allocate memory for each row */
    for (i = 0; i < n; i++) {
        A[i] = (double*)malloc(n * sizeof(double));
        if (A[i] == NULL) {
            free_matrix(A, i);
            exit_with_error();
        }
    }
    /* Compute similarity matrix values */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                dist = squared_euclidean_dist(datapoints[i], datapoints[j], d);
                A[i][j] = exp(-dist/2);
            } else {
                A[i][j] = 0.0;
            }
        }
    }
    return A;
}

/*
Given a pointer to a matrix M and the amount of rows in it, frees its memory (Every allocated row, and the row pointers array).
Assumes the pointer M leads to a positive-sized array of pointers.
*/
void free_matrix(double** M, int row_num)
{
    int i;
    for (i=0; i < row_num; i++)
    {
        if(M[i] != NULL)
            free(M[i]);
    }
    free(M);
}

/*
To be called in case of an error.
Prints an error message and terminates the program.
Note: FREE ALL DYNAMIC MEMORY BEFORE CALLING THIS FUNCTION!
*/
void exit_with_error() 
{
    printf(ERROR_MSG);
    exit(1);
}

/*
Given n*m matrix A, the value of m, a m*k matrix B, and two indices 0<=i<n, 0<=j<k.
Returns the value in cell (i,j) of the matrix AB.
This function is used to calculate H every iteration to save memory, instead of allocating and freeing lots of matrices every time.
*/
double matrix_mult_cell(double** A, int A_cols_num ,double** B, int i, int j) 
{   
    int p;
    double val = 0;
    for(p=0; p < A_cols_num; p++)
        val += (A[i][p] * B[p][j]);
    return val;
}

/*
Given two points represented as double-lists, return their Euclidean distance.
Assumes both points have the same dimension.
*/
double squared_euclidean_dist(double* point1, double* point2, int dimension)
{
    double sum = 0;
    int i;
    for (i = 0; i < dimension; i++)
    {
        sum += pow(point1[i] - point2[i], 2);
    }
    return sum;
}

/*
Given two NON-EMPTY matrices A,B and A's dimensions, calculates the Frobenius norm of A-B.
Assumes both matrices have the same dimensions.
*/
double frobenius_norm(double** A, int rows_num, int cols_num, double** B)
{ 
    int i,j;
    double sum = 0;
    for(i=0; i < rows_num; i++)
        for(j=0; j < cols_num; j++)
            sum += pow(A[i][j] - B[i][j], 2);
    return sqrt(sum);
}

/*
Given a n*m matrix M, its amount of rows and an int 0<=j<m, returns a 1*n matrix consisting only of M's j-th column.
If memory allocation error occurs, returns a null pointer.
*/
double** get_column(double** M, int rows_num, int j)
{   
    int i;
    double** ret = (double**)malloc(1 * sizeof(double*));
    if (ret == NULL)
        return NULL;
    ret[0] = (double*)malloc(rows_num * sizeof(double));
    if (ret[0] == NULL)
    {
        free_matrix(ret, 0);
        return NULL;
    }
    for(i=0; i < rows_num; i++)
        ret[0][i] = M[i][j];
    return ret;

}

/*
Given a n*n graph laplacian W, a current n*k iteration matrix H and a pointer to an ALREADY EXISTING n*k matrix new_H,
changes the values in the new_H matrix IN PLACE to be the new values, as per the instructions
(See 1.4.2).
If memory allocation error occurs, returns 1. if finished successfully, returns 0.
*/
int update_H(double** W, double** H, double** new_H, int n, int k)
{
    double numerator, denominator, cell_multiplier;
    double** Ht_row; /* The needed row in H^T to calculate the matrix product (H^T)H. Will be a 1*n matrix.*/
    int i,j,s;
    double** HtH_col = (double**)malloc(k * sizeof(double*)); /* The needed column in (H^T)H to calculate the denominator. Is a k*1 matrix. */
    if(HtH_col == NULL)
        return 1;
    for(j=0; j<k; j++)
    {
        HtH_col[j] = (double*)malloc(sizeof(double));
        if(HtH_col[j] == NULL)
        {
            free_matrix(HtH_col, j);
            return 1;
        }
    }
    for(j=0; j<k; j++)
    {
        for(s=0; s<k; s++) /* Calculates the necessary column of (H^T)H for the denominator. */
        {
            Ht_row = get_column(H, n, s);
            if(Ht_row == NULL)
            {
                free_matrix(HtH_col, k);
                return 1;
            }
            HtH_col[s][0] = matrix_mult_cell(Ht_row, n, H, 0, j);
            free_matrix(Ht_row, 1);
        }
        for(i=0; i<n; i++)
        {
            numerator = matrix_mult_cell(W, n, H, i, j);
            denominator = matrix_mult_cell(H, k, HtH_col, i, 0);
            cell_multiplier = numerator / denominator;
            cell_multiplier += (1 - beta);
            new_H[i][j] = H[i][j]*cell_multiplier;
        }
    }
    return 0;
}

/*
Given a starting matrix H, its dimensions and a graph laplacian W, perform the optimization algorithm in the instructions.
Returns an optimized H.
*/
double** optimizing_H(double** H, int rows_num, int cols_num, double** W)
{
    int i, j;
    double** tmp;
    double** new_H = (double**)malloc(rows_num * sizeof(double*));
    if (new_H == NULL)
    {
        free_matrix(H, rows_num);
        exit_with_error();
    }
    for (j=0; j < rows_num; j++)
    {
        new_H[j] = (double*)malloc(cols_num * sizeof(double));
        if(new_H[j] == NULL)
        {
            free_matrix(H, rows_num);
            free_matrix(new_H, j);
            exit_with_error();
        }
    }
    /* Does the actual work */
    for (i=1; i<=max_iter; i++)
    {
        if(update_H(W, H, new_H, rows_num, cols_num) == 1) /* Updates H and puts the updated version into new_H. 1 will be returned iff an error occurs during the update. */
        {
            free_matrix(H, rows_num);
            free_matrix(new_H, rows_num);
            exit_with_error();
        }
        if(frobenius_norm(new_H, rows_num, cols_num , H) < eps) /* We have reached convergence - end the loop. */
            i = max_iter + 1;
        tmp = H; /* Always makes the new matrix be in pointer H for code consistency. */
        H = new_H;
        new_H = tmp;
    }
    free_matrix(new_H, rows_num);
    return H;
}

/*
Receives a m*n matrix A and a n*k matrix B alongside their dimensions, and returns the product matrix AB.
*/
double** multiplyMatrix(double** matrixA, double** matrixB, int m, int n, int k){
    double** product = calloc(m, sizeof(double*));
    int i, j, l;
    for(i = 0; i < m; i++){
        product[i] = calloc(k, sizeof(double));
    }

    for (i = 0; i < m; i++){
        for (j = 0; j < k; j++){
            for (l = 0; l < n; l++){
                product[i][j] += matrixA[i][l] * matrixB[l][j];
            }
        }
    }

    return product;
}

/*
Receives two points (Represented as vectors) and their size, and returns the square of the euclidian norm of the vectors' difference.
Assumes both vectors are of the same dimension.
*/
double sq_euclidian_norm(double* point1, double* point2, int dimension){
    double sum = 0;
    int i;
    for (i = 0; i < dimension; i++){
        sum += pow(*(point1 + i) - *(point2 + i), 2);
    }
    return sum;
}

/*
Given an array of arrays representing points, the amount of points (n) and the dimension of every point (d),
returns the n*n similarity matrix of the points.
Assumes all points are of dimension d.
*/
double** similarity_matrix(double** datapoints, int n, int d){
    double** A = malloc(n * sizeof(double*));
    int i, j;
    if(A == NULL)
        exit_with_error();
    for (i = 0; i < n; i++){
        A[i] = calloc(n, sizeof(double));
        if(A[i] == NULL)
        {
            free_matrix(A, i);
            exit_with_error();
        }
    }
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i != j){
                A[i][j] = exp(-sq_euclidian_norm(*(datapoints + i), *(datapoints + j), d) / 2);
            } else {
                A[i][j] = 0;
            }
        }
    }
    return A;
}

/*
Given an n*n similarity matrix and the value of n, returns the normalized similarity matrix.
*/
double** normalized_similarity_matrix(double** sim_matrix, int n){
    int i;
    double** D = diagonal_degree_matrix(sim_matrix, n);
    double** D_neg_half = malloc(n * sizeof(double*));
    if(D_neg_half == NULL)
    {
        free_matrix(D, n);
        exit_with_error();
    }
    for (i = 0; i < n; i++){
        D_neg_half[i] = calloc(n, sizeof(double));
        if(D_neg_half[i] == NULL)
        {
            free_matrix(D, n); free_matrix(D_neg_half, i);
            exit_with_error();
        }
    }
    for (i = 0; i < n; i++){
        D_neg_half[i][i] = 1 / sqrt(D[i][i]);
    }
    double** temp = multiplyMatrix(D_neg_half, sim_matrix, n, n, n);
    double** normalized = multiplyMatrix(temp, D_neg_half, n, n, n);
    free_matrix(D, n);
    free_matrix(D_neg_half, n);
    free_matrix(temp, n);
    return normalized;
}