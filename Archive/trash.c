#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

/* OLD initH in py
# def initH(n, k, W):
#     print(n, k)
#     if isinstance(W, list):
#         W = pd.DataFrame(W)
#         m = W.values.mean()
#         H = pd.DataFrame(np.random.uniform(0, 2*math.sqrt(m/k), size=(n, k)))
#         return H
*/