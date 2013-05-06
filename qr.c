#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
    Modified Gram-Schmidt algorithm as given by Trefethen
    Row-major to optimize cache access
*/

double
innerproduct_dd(double *v1, double *v2, size_t m)
{
    double result = 0.0;
    int i;
    for (i = 0; i < m; ++i)
    {
        result += v1[i] * v2[i];
    }
    return result;
}

void
gram_schmidt_d(double **rop, double **input, size_t m, size_t n)
{
    int i, j, k;
    for (i = 0; i < n; ++i)
    {
        double norm;
        norm = 0.0;
        for (j = 0; j < m; ++j)
        {
            rop[i][j] = input[i][j];
            norm += rop[i][j] * rop[i][j];
        }
        norm = sqrt(norm);
        if (norm > 1e-13)
        {
            for (j = 0; j < m; ++j)
            {
                rop[i][j] = rop[i][j] / norm;
            }
        }
        for (k = i + 1; k < n; ++k)
        {
            norm = innerproduct_dd(rop[i], rop[k], m);
            for (j = 0; j < m; ++j)
            {
                rop[k][j] = input[k][j] - norm * rop[i][j];
            }
        }
    }
}

void
print_matrix_d(double **matrix, size_t m, size_t n)
{
    int i, j;
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            printf("%.5f ", matrix[i][j]);
        }
        printf("\n");
    }
}

double **
alloc_matrix_d(size_t m, size_t n)
{
    int i, j;
    double **A;
    A = calloc(m, sizeof(double *));
    for (i = 0; i < m; ++i)
    {
        A[i] = calloc(n, sizeof(double));
        for (j = 0; j < n; ++j)
        {
            scanf("%20lf", &A[i][j]);
        }
    }
    return A;
}

void
free_matrix_d(double **A, size_t m)
{
    int i;
    for (i = 0; i < m; ++i)
    {
        free(A[i]);
    }
    free(A);
}

int
main(void)
{
    int m, n;
    double **A;
    scanf("%20d", &m);
    scanf("%20d", &n);

    if (m < 1 || n < 1)
    {
        printf("The input matrix should be at least 1x1\n");
        return -1;
    }

    A = alloc_matrix_d(m, n);

    gram_schmidt_d(A, A, m, n);
    print_matrix_d(A, m, n);

    free_matrix_d(A, m);

    return 0;
}
