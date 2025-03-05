/*
* symnmf.c - C interface for symNMF implementation
* This file provides an executable interface for the symNMF functions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "symnmf.h"

#define SEPARATOR ","
#define ERROR_MSG "An Error Has Occurred\n"
#define MAX_LINE_LENGTH 1024

/* Function declarations */
double **read_data(const char *filename, int *n, int *d);
void free_matrix(double **data, int n);
void print_matrix(double **matrix, int rows, int cols);

/* 
read_data - Reads data points from a file
Parameters:
    filename - Path to the input file
    n - Pointer in which to store the number of data points
    d - Pointer in which to store the dimension of each data point
Returns: Double pointer to the data points matrix (n x d)
*/
double **read_data(const char *filename, int *n, int *d) {
    FILE *fp;
    double **points;
    char line[MAX_LINE_LENGTH];
    char *token;
    int i, j;
    
    /* Open file */
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s\n", ERROR_MSG);
        exit(1);
    }
    
    *n = 0; *d = 0; /* Count num of points and dimensions */
    
    /* Read first line to count dimensions (=d) */
    if (fgets(line, sizeof(line), fp) != NULL) {
        token = strtok(line, SEPARATOR); /* Like "split" in py */
        while (token != NULL) {
            (*d)++;
            token = strtok(NULL, SEPARATOR);
        }
        (*n)++;
    }
    /* Then count remaining lines (=points=n) */
    while (fgets(line, sizeof(line), fp) != NULL) {
        (*n)++;
    }
    
    /* Reset file pointer to beginning of file */
    rewind(fp);
    
    /* Allocate memory for data points matrix */
    points = (double **)malloc(*n * sizeof(double *));
    if (points == NULL) {
        printf("%s\n", ERROR_MSG);
        fclose(fp);
        exit(1);
    }
    
    /* Allocate memory for each row */
    for (i = 0; i < *n; i++) {
        points[i] = (double *)malloc(*d * sizeof(double));
        if (points[i] == NULL) {
            printf("%s\n", ERROR_MSG);
            /* If theres a problem in allocation, 
            free previously allocated memory and exit */
            for (j = 0; j < i; j++) {
                free(points[j]);
            }
            free(points);
            fclose(fp);
            exit(1);
        }
    }
    
    /* Read data points from file */
    i = 0;
    while (fgets(line, sizeof(line), fp) != NULL && i < *n) {
        token = strtok(line, SEPARATOR); /* Like "split" in py */
        j = 0;
        while (token != NULL && j < *d) {
            points[i][j] = atof(token);
            token = strtok(NULL, SEPARATOR); /* String to float */
            j++;
        }
        i++;
    }
    
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
 * CMD args:
 * argv[1] - goal (sym, ddg, or norm)
 * argv[2] - file path
 */
int main(int argc, char *argv[]) {
    double **points, **result = NULL;
    int n, d;
    char *goal = argv[1];
    char *filename = argv[2];
    
    /* Check for correct numb of CMD args */
    if (argc != 3) {
        printf("%s %s %s\n", argv[0], argv[1], argv[2]);
        printf("%s\n", ERROR_MSG);
        return 1;
    }
    
    /* Read data points from input file */
    points = read_data(filename, &n, &d);
    
    // print_matrix(points, n, d);
    // exit(0);

    /* Call appropriate function based on goal */
    // if (strcmp(goal, "sym") == 0) {
    //     /* Calculate similarity matrix */
    //     result = sym(points, n, d);
    //     if (result == NULL) {
    //         printf("%s\n", ERROR_MSG);
    //         free_matrix(points, n);
    //         return 1;
    //     }
    //     print_matrix(result, n, n);
    // } else if (strcmp(goal, "ddg") == 0) {
    //     /* Calculate diagonal degree matrix */
    //     result = ddg(points, n, d);
    //     if (result == NULL) {
    //         printf("%s\n", ERROR_MSG);
    //         free_matrix(points, n);
    //         return 1;
    //     }
    //     print_matrix(result, n, n);
    // } else if (strcmp(goal, "norm") == 0) {
    //     /* Calculate normalized similarity matrix */
    //     result = norm(points, n, d);
    //     if (result == NULL) {
    //         printf("%s\n", ERROR_MSG);
    //         free_matrix(points, n);
    //         return 1;
    //     }
    //     print_matrix(result, n, n);
    // } else {
    //     /* Invalid goal */
    //     printf("%s\n", ERROR_MSG);
    //     free_matrix(points, n);
    //     return 1;
    // }
    
    /* free allocated memory */
    free_matrix(points, n);
    free_matrix(result, n);
    return 0;
}