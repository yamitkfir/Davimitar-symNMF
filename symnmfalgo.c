#include <math.h>

#define max_iter 300
#define eps 1e-4
#define beta 0.5
/* the actual mess */

/* 
similarity_matrix()
diagonal_matrix()
normalized_similarity_matrix()
init_H()
*/

/* Function declarations */
void free_matrix();
void exit_with_error();
double matrix_mult_cell();
double euclidian_dist();
double frobenius_norm();
double** get_column();
double** optimizing_H();
int update_H();

void free_matrix(double** M)
/*
Given a pointer to a matrix M, frees its memory.
Assumes the pointer M leads to a positive-sized array of pointers.
*/
{
    int i;
    int len = sizeof(M) / sizeof(double*);
    for (i=0; i<len; i++)
        if(M[i] != 0)
            free(M[i]);
    free(M);
}

void exit_with_error() // Developer's note: FREE ALL DYNAMIC MEMORY BEFORE CALLING THIS FUNCTION!!!
/*
To be called in case of an error.
Prints an error message and terminates the program.
*/
{
    printf("An error has occurred.");
    exit();
}

double matrix_mult_cell(double** A, double** B, int i, int j)
/*
Given n*m matrix A, m*k matrix B, and two indices 0<=i<n, 0<=j<k.
Returns the value in cell (i,j) of the matrix AB.
*/
{
    int A_cols_num = sizeof(A[0]) / sizeof(double);
    int p;
    double val = 0;
    for(p=0; p<A_cols_num; p++)
        val += (A[i][p] * B[p][j]);
    return val;
}

double euclidian_dist(){
    /* TODO yamit - We have to agree on a format for points before this one! */
    return 0;
}

double frobenius_norm(double** A, double** B)
/*
Given two NON-EMPTY matrices A,B, calculates the Frobenius norm of A-B.
Assumes both matrices have the same dimensions.
*/
{
    int rows_num = sizeof(A) / sizeof(double*);
    int cols_num = sizeof(A[0]) / sizeof(double);
    int i,j;
    double sum = 0;
    for(i=0; i<rows_num; i++)
        for(j=0; j<cols_num; j++)
            sum += pow(A[i][j]-B[i][j], 2);
    return sqrt(sum);
}

double** get_column(double** M, int j)
/*
Given a n*m matrix M and an int 0<=j<m, returns a 1*n matrix consisting only of M's j-th column.
If memory allocation error occurs, returns a null pointer.
*/
{
    int rows_num = sizeof(M) / sizeof(double*);
    int i;
    double** ret = (double**)calloc(sizeof(double*));
    if (ret == 0)
        return 0;
    ret[0] = (double*)malloc(rows_num * sizeof(double));
    if (ret[0] == 0)
    {
        free_matrix(ret);
        return 0;
    }
    for(i=0; i<rows_num; i++)
        ret[0][i] = M[i][j];
    return ret;

}

int update_H(double** W, double** H, double** new_H, int n, int k)
/*
Given a n*n graph laplacian W, a current n*k iteration matrix H and a pointer to an ALREADY EXISTING n*k matrix new_H,
changes the values in the new_H matrix IN PLACE to be the new values, as per the instructions
(See 1.4.2).
If memory allocation error occurs, returns 1. if finished successfully, returns 0.
*/
{
    double numerator, denominator, cell_multiplier;
    double** Ht_row; // The needed row in H^T to calculate the matrix product (H^T)H.
    int i,j,s;
    // The needed column in (H^T)H to calculate the denominator. Is a k*1 matrix.
    double** HtH_col = (double**)calloc(k*sizeof(double*));
    if(HtH_col == 0)
        return 1;
    for(j=0; j<k; j++)
    {
        HtH_col[j] = (double*)malloc(sizeof(double));
        if(HtH_col[0][j] == 0)
        {
            free_matrix(HtH_col);
            return 1;
        }
    }
    for(j=0; j<k; j++)
    {
        for(s=0; s<k; s++) // Calculates the necessary column of (H^T)H for the denominator.
        {
            Ht_row = get_column(H, s);
            if(Ht_row == 0)
            {
                free_matrix(HtH_col);
                return 1;
            }
            HtH_col[s][0] = matrix_mult_cell(Ht_row, H, 0, j);
            free_matrix(Ht_row);
        }
        for(i=0; i<n; i++)
        {
            numerator = matrix_mult_cell(W, H, i, j);
            denominator = matrix_mult_cell(H, HtH_col, i, 0);
            cell_multiplier = numerator / denominator;
            cell_multiplier += (1 - beta);
            new_H[i][j] = H[i][j]*cell_multiplier;
        }
    }
    return 0;
}

double** optimizing_H(double** H, double** W)
/*
Given a starting matrix H and a graph laplacian W, perform the optimization algorithm in the instructions.
Returns an optimized H.
*/
{
    int rows_num = sizeof(H) / sizeof(double*); // Number of rows in H
    int cols_num = sizeof(H[0]) / sizeof(double); // Number of columns in H
    int i, j;
    double** tmp;
    double** new_H = (double**)calloc(rows_num * sizeof(double*));
    if (new_H == 0)
    {
        free_matrix(H);
        exit_with_error();
    }
    for (j=0; j<rows_num; i++)
    {
        new_H[j] = (double*)malloc(cols_num * sizeof(double));
        if(new_H[j] == 0)
        {
            free_matrix(H);
            free_matrix(new_H);
            exit_with_error();
        }
    }
    // Does the actual work
    for (i=1; i<=max_iter; i++)
    {
        if(update_H(W, H, new_H, rows_num, cols_num) == 1)
        {
            free_matrix(H);
            free_matrix(new_H);
            exit_with_error();
        }
        if(frobenius_norm(new_H, H) < eps) // We have reached convergence - end the loop.
            i = max_iter + 1;
        tmp = H; // Always makes the new matrix be in pointer H for code consistency.
        H = new_H;
        new_H = tmp;
    }
    free_matrix(new_H);
    return H;
}