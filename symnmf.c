/*
* symnmf.c - C interface for symNMF implementation
* This file provides an executable interface for the symNMF functions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* #include "symnmf.h" */

#define max_iter 300
#define eps 1e-4
#define denominator_eps 1e-7
#define beta 0.5
#define SEPARATOR ","
#define ERROR_MSG "An Error Has Occurred\n"
#define MAX_LINE_LENGTH 1024

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
double **create_points_matrix(FILE *fp, char line[], int *n, int *d);
double frobenius_norm(double** A, int rows_num, int cols_num, double** B);
void free_matrix(double** M, int len);
double** multiply_matrix(double** matrixA, double** matrixB, int m, int n, int k); /* A - m x n, B - n x k */
double matrix_mult_cell(double** A, int A_cols_num ,double** B, int i, int j);
void update_H_cell(double **W, double **H, double **new_H, double **HtH_col, int n, int k, int i, int j);
void exit_with_error();


void exit_with_error() 
/* note: FREE ALL DYNAMIC MEMORY BEFORE CALLING THIS FUNCTION! */
{
    printf(ERROR_MSG);
    exit(1);
}

/*
Given an opened file fp, an array big enough to hold every line from fp and the dimensions of the points represented in fp, returns a n*d point matrix of the points in the file.
*/
double **create_points_matrix(FILE *fp, char line[], int *n, int *d)
{
    int i, j;
    char *token;
    double **points = (double **)malloc(*n * sizeof(double *)); /* Allocate memory for data points matrix */
    if (points == NULL) {
        fclose(fp);
        exit_with_error();
    }
    for (i = 0; i < *n; i++) { /* Allocate memory for each row */
        points[i] = (double *)malloc(*d * sizeof(double));
        if (points[i] == NULL) { /* If theres a problem in allocation, free previously allocated memory and exit */
            for (j = 0; j < i; j++) 
                free(points[j]);
            free(points);
            fclose(fp);
            exit_with_error();
        }
    }
    i = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL && i < *n) { /* Read data points from file */
        token = strtok(line, SEPARATOR); /* Like "split" in Python */
        j = 0;
        while (token != NULL && j < *d) {
            points[i][j] = atof(token);
            token = strtok(NULL, SEPARATOR); /* String to float */
            j++;
        }
        i++;
    }
    return points;
}

/* 
Reads data points from a file.
Parameters: filename - Path to the input file, n - Pointer in which to store the number of data points, d - Pointer in which to store the dimension of each data point
Returns: Double pointer to the data points matrix (n x d)
*/
double **read_data(const char *filename, int *n, int *d) {
    FILE *fp;
    double **points;
    char line[MAX_LINE_LENGTH];
    char *token;
    fp = fopen(filename, "r"); /* Open file */
    if (fp == NULL) {
        exit_with_error();
    }
    *n = 0; *d = 0; /* Count num of points and dimensions */
    if (fgets(line, MAX_LINE_LENGTH, fp) != NULL) { /* Read first line to count dimensions (=d) */
        token = strtok(line, SEPARATOR); /* Like "split" in py */
        while (token != NULL) {
            (*d)++;
            token = strtok(NULL, SEPARATOR);
        }
        (*n)++;
    }
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) /* Then count remaining lines (=points=n) */
        (*n)++;
    rewind(fp); /* Reset file pointer to beginning of file */
    points = create_points_matrix(fp, line, n, d);
    fclose(fp);
    return points;
}

/*
* free_matrix - Frees memory allocated for a matrix
* n - Number of rows
*/
void free_matrix(double **data, int n) {
    int i;
    if (data == NULL) {
        return;
    }
    for (i = 0; i < n; i++) {
        if (data[i] != NULL) {
            free(data[i]);
        }
    }
    free(data);
}

void print_matrix(double **matrix, int rows, int cols) {
    int i, j;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < cols - 1) {
                printf("%s", SEPARATOR);
            }
        }
        printf("\n");
    }
}


/*
Calculate the Diagonal Degree Matrix D for a given similarity matrix A.
Uses a 2D array representation (array of arrays).
*/
double** diagonal_degree_matrix(double** A, int n) {
    double** D = (double**)malloc(n * sizeof(double*));
    int i, j; double sum;
    if (D == NULL) {
        exit_with_error();
    }
    /* Allocate memory for each row and init to 0 */
    for (i = 0; i < n; i++) {
        D[i] = (double*)malloc(n * sizeof(double)); /* TODO david replaced calloc with malloc, ok? */
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
Given the current iteration's H matrix, a pointer to the new H matrix and all needed values (including a cell coordinate (i,j) in H), updates the new H's cell in place (i,j).
*/
void update_H_cell(double **W, double **H, double **new_H, double **HtH_col, int n, int k, int i, int j)
{
    int numerator, denominator, cell_multiplier;
    numerator = matrix_mult_cell(W, n, H, i, j);
    denominator = matrix_mult_cell(H, k, HtH_col, i, 0) + denominator_eps; /* This epsilon is added to avoid division by zero. */
    cell_multiplier = numerator / denominator;
    cell_multiplier += (1 - beta);
    new_H[i][j] = H[i][j]*cell_multiplier;
}

/*
Given a n*n graph laplacian W, a current n*k iteration matrix H and a pointer to an ALREADY EXISTING n*k matrix new_H,
changes the values in the new_H matrix IN PLACE to be the new values, as per the instructions (See 1.4.2).
If memory allocation error occurs, returns 1. if finished successfully, returns 0.
*/
int update_H(double** W, double** H, double** new_H, int n, int k)
{
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
            update_H_cell(W, H, new_H, HtH_col, n, k, i, j);
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
    double **tmp, **new_H = (double**)malloc(rows_num * sizeof(double*));
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
    for (i=1; i<=max_iter; i++) /* Does the actual work */
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
double** multiply_matrix(double** matrixA, double** matrixB, int m, int n, int k) {
    double** product = malloc(m * sizeof(double*));
    if (product == NULL) return NULL;
    int i, j, l;
    for (i = 0; i < m; i++) {
        product[i] = malloc(k * sizeof(double));
        if (product[i] == NULL) {
            // Clean up previously allocated memory
            for (int p = 0; p < i; p++) {
                free(product[p]);
            }
            free(product);
            return NULL;
        }
        /* Initialize to zero (still more efficient than calloc) */
        for (j = 0; j < k; j++) {
            product[i][j] = 0.0;
        }
    }
    /* Optimized loop order for better cache locality */
    for (i = 0; i < m; i++) {
        for (l = 0; l < n; l++) {
            const double a_il = matrixA[i][l]; /* Hear me out: Cache this value */
            for (j = 0; j < k; j++) {
                product[i][j] += a_il * matrixB[l][j];
            }
        }
    }
    return product;
}

/*
Given an array of arrays representing points, the amount of points (n) and the dimension of every point (d),
returns the n*n similarity matrix of the points. Assumes all points are of dimension d.
TODO free_matrix(A) once done with it!
*/
double** similarity_matrix(double** datapoints, int n, int d){
    double** A = malloc(n * sizeof(double*));
    int i, j;
    if(A == NULL)
        return NULL;
    for (i = 0; i < n; i++){
        A[i] = malloc(n * sizeof(double)); /* TODO david replaced calloc with malloc, ok? */
        if(A[i] == NULL)
        {
            free_matrix(A, i);
            return NULL;
        }
    }
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            if (i != j){
                A[i][j] = exp(-squared_euclidean_dist(*(datapoints + i), *(datapoints + j), d) / 2);
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
    double** temp, **normalized;
    double** D = diagonal_degree_matrix(sim_matrix, n);
    double** D_neg_half = malloc(n * sizeof(double*));
    if(D_neg_half == NULL)
    {
        free_matrix(D, n);
        exit_with_error();
    }
    for (i = 0; i < n; i++){
        D_neg_half[i] = malloc(n * sizeof(double)); /* TODO david replaced calloc with malloc, ok? */
        if(D_neg_half[i] == NULL)
        {
            free_matrix(D, n); free_matrix(D_neg_half, i);
            exit_with_error();
        }
    }
    for (i = 0; i < n; i++){
        D_neg_half[i][i] = 1 / sqrt(D[i][i]);
    }
    printf("425\n"); /* TEMP */
    temp = multiply_matrix(D_neg_half, sim_matrix, n, n, n);
    printf("427\n"); /* TEMP */
    normalized = multiply_matrix(temp, D_neg_half, n, n, n);
    printf("429\n"); /* TEMP */
    free_matrix(D, n);
    free_matrix(D_neg_half, n);
    free_matrix(temp, n);
    return normalized;
}


/*
CMD args: argv[1] - goal (sym, ddg, or norm), argv[2] - file path
*/
int main(int argc, char *argv[]) {
    double **points, **A = NULL , **result = NULL;
    int n, d;
    char *goal = argv[1], *filename = argv[2];
    if (argc != 3) { exit_with_error(); } /* Check for correct num of CMD args */
    points = read_data(filename, &n, &d); /* Read data points from input file */
    if (strcmp(goal, "sym") == 0) { /* Goal: calculate similarity matrix */
        result = similarity_matrix(points, n, d); 
        if (result == NULL) {
            free_matrix(points, n);
            exit_with_error();
        }
    } else if (strcmp(goal, "ddg") == 0) { /* Goal: calculate diagonal degree matrix */
        A = similarity_matrix(points, n, d);
        if (A == NULL) { /* Couldn't allocate space for A */
            free_matrix(points, n);
            exit_with_error();
        }
        result = diagonal_degree_matrix(A, n);
        if (result == NULL) {
            free_matrix(points, n);
            free_matrix(A, n);
            exit_with_error();
        }
    } else if (strcmp(goal, "norm") == 0) { /* Goal: calculate normalized similarity matrix */
        A = similarity_matrix(points, n, d);
        if (A == NULL) { /* Couldn't allocate space for A */
            free_matrix(points, n);
            exit_with_error();
        }
        result = normalized_similarity_matrix(A, n);
        if (result == NULL) {
            free_matrix(points, n);
            free_matrix(A, n);
            exit_with_error();
        }  
    } else { /* Invalid goal */
        free_matrix(points, n);
        exit_with_error();
    }
    print_matrix(result, n, n);
    free_matrix(points, n);
    free_matrix(result, n);
    if (strcmp(goal, "sym") != 0) { /* Only free A if it was actually allocated */
        free_matrix(A, n);
    }
    return 0;
}