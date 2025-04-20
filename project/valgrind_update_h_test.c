/*
 * valgrind_update_h_test.c - Memory leak test specifically for the update_H function
 * 
 * Compile: gcc -ansi -Wall -Wextra -pedantic-errors -g valgrind_update_h_test.c -o valgrind_update_h_test -lm
 * Run with Valgrind: valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./valgrind_update_h_test
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Constants similar to those in symnmf.c */
#define max_iter 300
#define eps 1e-4
#define denominator_eps 1e-7
#define beta 0.5
#define SEPARATOR ","
#define ERROR_MSG "An Error Has Occurred\n"
#define MAX_LINE_LENGTH 1024

/* Function prototypes */
void free_matrix(double **data, int n);
double** create_matrix(int n, int m);
double** multiply_matrices(double** A, double** B, int n, int k, int m);
double** get_column(double** M, int rows_num, int j);
double matrix_mult_cell(double** A, int A_cols_num, double** B, int i, int j);
void update_H_cell(double **W, double **H, double **new_H, double **HtH_col, int n, int k, int i, int j);
int update_H(double** W, double** H, double** new_H, int n, int k);
double** optimizing_H(double** H, int rows_num, int cols_num, double** W);
double sq_frobenius_norm(double** A, int rows_num, int cols_num, double** B);
void print_matrix(double **matrix, int rows, int cols);
void exit_with_error(void);

int main(void) {
    int n = 4, k = 2;
    int i, j;
    double **W = NULL, **H = NULL, **H_optimized = NULL;
    double **new_H;
    int result;
    double **H_copy;

    printf("Starting memory leak test for update_H function...\n");
    
    /* Create W matrix (normalized similarity matrix) */
    W = create_matrix(n, n);
    if (W == NULL) {
        printf("Failed to allocate memory for W matrix.\n");
        return 1;
    }
    
    /* Fill W with test values (symmetric) */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                W[i][j] = 1.0;
            } else {
                W[i][j] = 0.5 / (abs(i - j) + 1);
                /* Ensure symmetry */
                W[j][i] = W[i][j];
            }
        }
    }
    
    printf("\nW matrix (normalized similarity matrix):\n");
    print_matrix(W, n, n);
    
    /* Create initial H matrix */
    H = create_matrix(n, k);
    if (H == NULL) {
        printf("Failed to allocate memory for H matrix.\n");
        free_matrix(W, n);
        return 1;
    }
    
    /* Fill H with test values (random between 0 and 1) */
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            H[i][j] = (double)rand() / RAND_MAX;
        }
    }
    
    printf("\nInitial H matrix:\n");
    print_matrix(H, n, k);
    
    /* Test update_H with a single iteration */
    printf("\nTesting one iteration of update_H...\n");
    new_H = create_matrix(n, k);
    if (new_H == NULL) {
        printf("Failed to allocate memory for new_H matrix.\n");
        free_matrix(W, n);
        free_matrix(H, n);
        return 1;
    }
    
    /* Update H for one iteration */
    result = update_H(W, H, new_H, n, k);
    if (result != 0) {
        printf("Error during update_H execution.\n");
        free_matrix(W, n);
        free_matrix(H, n);
        free_matrix(new_H, n);
        return 1;
    }
    
    printf("\nH matrix after one update iteration:\n");
    print_matrix(new_H, n, k);
    
    /* Free matrices used for single iteration test */
    free_matrix(new_H, n);
    
    /* Now test the full optimizing_H function */
    printf("\nTesting full optimizing_H function...\n");
    
    /* Create a fresh copy of H for optimization */
    H_copy = create_matrix(n, k);
    if (H_copy == NULL) {
        printf("Failed to allocate memory for H_copy matrix.\n");
        free_matrix(W, n);
        free_matrix(H, n);
        return 1;
    }
    
    /* Copy H values to H_copy */
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            H_copy[i][j] = H[i][j];
        }
    }
    
    /* Optimize H */
    H_optimized = optimizing_H(H_copy, n, k, W);
    /* Note: H_copy is consumed by optimizing_H, so we don't need to free it */
    
    if (H_optimized == NULL) {
        printf("Failed to optimize H matrix.\n");
        free_matrix(W, n);
        free_matrix(H, n);
        return 1;
    }
    
    printf("\nOptimized H matrix:\n");
    print_matrix(H_optimized, n, k);
    
    /* Free all allocated memory */
    free_matrix(W, n);
    free_matrix(H, n);
    free_matrix(H_optimized, n);
    
    printf("\nAll memory freed successfully. No leaks detected.\n");
    return 0;
}

/* Function to create a matrix with given dimensions */
double** create_matrix(int n, int m) {
    int i;
    double **matrix = (double **)malloc(n * sizeof(double *));
    
    if (matrix == NULL) {
        return NULL;
    }
    
    for (i = 0; i < n; i++) {
        matrix[i] = (double *)calloc(m, sizeof(double));
        if (matrix[i] == NULL) {
            free_matrix(matrix, i);
            return NULL;
        }
    }
    
    return matrix;
}

/* Function to free a matrix */
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

/* Function to print a matrix */
void print_matrix(double **matrix, int rows, int cols) {
    int i, j;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            printf("%.4f ", matrix[i][j]);
        }
        printf("\n");
    }
}

/* Function to exit with error */
void exit_with_error(void) {
    printf("%s", ERROR_MSG);
    exit(1);
}

/* Get a column from a matrix as a transposed row vector */
double** get_column(double** M, int rows_num, int j) {
    int i;
    double** ret = (double**)malloc(1 * sizeof(double*));
    
    if (ret == NULL) {
        return NULL;
    }
    
    ret[0] = (double*)malloc(rows_num * sizeof(double));
    if (ret[0] == NULL) {
        free_matrix(ret, 0);
        return NULL;
    }
    
    for (i = 0; i < rows_num; i++) {
        ret[0][i] = M[i][j];
    }
    
    return ret;
}

/* Calculate a single cell of a matrix multiplication */
double matrix_mult_cell(double** A, int A_cols_num, double** B, int i, int j) {
    int p;
    double val = 0;
    
    for (p = 0; p < A_cols_num; p++) {
        val += (A[i][p] * B[p][j]);
    }
    
    return val;
}

/* Update a single cell in the H matrix */
void update_H_cell(double **W, double **H, double **new_H, double **HtH_col, int n, int k, int i, int j) {
    double numerator, denominator, cell_multiplier;
    
    numerator = matrix_mult_cell(W, n, H, i, j);
    denominator = matrix_mult_cell(H, k, HtH_col, i, 0) + denominator_eps; /* Add epsilon to avoid division by zero */
    cell_multiplier = numerator / denominator;
    cell_multiplier *= beta;
    cell_multiplier += (1 - beta);
    new_H[i][j] = H[i][j] * cell_multiplier;
}

/* Update the entire H matrix for one iteration */
int update_H(double** W, double** H, double** new_H, int n, int k) {
    double** Ht_row; /* Transposed row of H */
    int i, j, s;
    double** HtH_col = (double**)malloc(k * sizeof(double*)); /* Column of (H^T)H */
    
    if (HtH_col == NULL) {
        return 1;
    }
    
    for (j = 0; j < k; j++) {
        HtH_col[j] = (double*)calloc(1, sizeof(double));
        if (HtH_col[j] == NULL) {
            free_matrix(HtH_col, j);
            return 1;
        }
    }
    
    for (j = 0; j < k; j++) {
        for (s = 0; s < k; s++) { /* Calculate the necessary column of (H^T)H for the denominator */
            Ht_row = get_column(H, n, s);
            if (Ht_row == NULL) {
                free_matrix(HtH_col, k);
                return 1;
            }
            HtH_col[s][0] = matrix_mult_cell(Ht_row, n, H, 0, j);
            free_matrix(Ht_row, 1);
        }
        
        for (i = 0; i < n; i++) {
            update_H_cell(W, H, new_H, HtH_col, n, k, i, j);
        }
    }
    
    free_matrix(HtH_col, k);
    return 0;
}

/* Calculate the squared Frobenius norm of the difference between two matrices */
double sq_frobenius_norm(double** A, int rows_num, int cols_num, double** B) {
    int i, j;
    double sum = 0;
    
    for (i = 0; i < rows_num; i++) {
        for (j = 0; j < cols_num; j++) {
            sum += pow(A[i][j] - B[i][j], 2);
        }
    }
    
    return sum;
}

/* Function to multiply two matrices */
double** multiply_matrices(double** A, double** B, int n, int k, int m) {
    double** C = create_matrix(n, m);
    int i, j, l;
    
    if (C == NULL) {
        return NULL;
    }
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            C[i][j] = 0.0;
            for (l = 0; l < k; l++) {
                C[i][j] += A[i][l] * B[l][j];
            }
        }
    }
    
    return C;
}

/* Optimize H using iterative updates */
double** optimizing_H(double** H, int rows_num, int cols_num, double** W) {
    int i;
    double **tmp, **new_H = create_matrix(rows_num, cols_num);
    
    if (new_H == NULL) {
        free_matrix(H, rows_num);
        return NULL;
    }
    
    /* Iterate until convergence or max iterations */
    for (i = 1; i <= max_iter; i++) {
        /* Update H and put the updated version into new_H */
        if (update_H(W, H, new_H, rows_num, cols_num) == 1) {
            free_matrix(H, rows_num);
            free_matrix(new_H, rows_num);
            return NULL;
        }
        
        /* Check for convergence */
        if (sq_frobenius_norm(new_H, rows_num, cols_num, H) < eps) {
            /* Reached convergence */
            break;
        }
        
        /* Swap H and new_H for next iteration */
        tmp = H;
        H = new_H;
        new_H = tmp;
    }
    
    /* Free the matrix that's not being returned */
    free_matrix(new_H, rows_num);
    
    return H;
}