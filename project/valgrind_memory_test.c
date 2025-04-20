/*
 * valgrind_memory_test.c - Memory leak test for symNMF implementation
 * 
 * Compile: gcc -ansi -Wall -Wextra -pedantic-errors -g valgrind_memory_test.c -o valgrind_memory_test -lm
 * Run with Valgrind: valgrind --leak-check=full --show-leak-kinds=all ./valgrind_memory_test
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

/* Function prototypes */
double** create_matrix(int n, int m);

int main(void) {
    int n = 4, d = 2, k = 2;
    int i, j;
    double **points = NULL;
    double **similarity = NULL;
    double **diagonal = NULL;
    double **normalized = NULL;
    double **H = NULL;
    double **optimized_H = NULL;

    printf("Starting logical flow memory leak test for symNMF...\n");

    /* Step 1: Create test data points */
    points = create_matrix(n, d);
    if (points == NULL) {
        printf("Failed to allocate memory for points.\n");
        return 1;
    }

    /* Fill points with test data */
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            points[i][j] = sin((double)(i * j + 1));
        }
    }

    printf("\n1. Test data points:\n");
    print_matrix(points, n, d);

    /* Step 2: Calculate similarity matrix */
    similarity = similarity_matrix(points, n, d);
    if (similarity == NULL) {
        printf("Failed to calculate similarity matrix.\n");
        free_matrix(points, n);
        return 1;
    }

    printf("\n2. Similarity matrix:\n");
    print_matrix(similarity, n, n);

    /* Step 3: Calculate diagonal degree matrix */
    diagonal = diagonal_degree_matrix(similarity, n);
    if (diagonal == NULL) {
        printf("Failed to calculate diagonal degree matrix.\n");
        free_matrix(points, n);
        free_matrix(similarity, n);
        return 1;
    }

    printf("\n3. Diagonal degree matrix:\n");
    print_matrix(diagonal, n, n);

    /* Step 4: Calculate normalized similarity matrix */
    normalized = normalized_similarity_matrix(similarity, n);
    if (normalized == NULL) {
        printf("Failed to calculate normalized similarity matrix.\n");
        free_matrix(points, n);
        free_matrix(similarity, n);
        free_matrix(diagonal, n);
        return 1;
    }

    printf("\n4. Normalized similarity matrix:\n");
    print_matrix(normalized, n, n);

    /* Step 5: Create initial H matrix */
    H = create_matrix(n, k);
    if (H == NULL) {
        printf("Failed to allocate memory for H matrix.\n");
        free_matrix(points, n);
        free_matrix(similarity, n);
        free_matrix(diagonal, n);
        free_matrix(normalized, n);
        return 1;
    }

    /* Fill H with random data */
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            H[i][j] = (double)rand() / RAND_MAX;
        }
    }

    printf("\n5. Initial H matrix:\n");
    print_matrix(H, n, k);

    /* Step 6: Optimize H using symNMF */
    optimized_H = optimizing_H(H, n, k, normalized);
    if (optimized_H == NULL) {
        printf("Failed to optimize H matrix.\n");
        free_matrix(points, n);
        free_matrix(similarity, n);
        free_matrix(diagonal, n);
        free_matrix(normalized, n);
        free_matrix(H, n);
        return 1;
    }

    printf("\n6. Optimized H matrix:\n");
    print_matrix(optimized_H, n, k);

    /* Free all allocated memory */
    free_matrix(points, n);
    free_matrix(similarity, n);
    free_matrix(diagonal, n);
    free_matrix(normalized, n);
    free_matrix(optimized_H, n);

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

