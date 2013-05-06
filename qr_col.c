#include "qr.h"

/*
    Modified Gram-Schmidt algorithm as given by Golub
    Optimized to better use the cache
*/

double **
alloc_matrix_d(size_t m, size_t n)
{
    int j;
    double **A;
    A = calloc(n, sizeof(double *));
    for (j = 0; j < n; ++j)
    {
        A[j] = calloc(m, sizeof(double));
    }
    return A;
}

double
col_product_dd(double **A, double **B, size_t a, size_t b, size_t m)
{
    double result = 0.0;
    int i;
    for (i = 0; i < m; ++i)
    {
        result += A[a][i] * B[b][i];
    }
    return result;
}

void
copy_matrix_dd(double **rop, double **A, size_t m, size_t n)
{
    int i, j;
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            rop[j][i] = A[j][i];
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
            printf("%10.5f ", matrix[j][i]);
        }
        printf("\n");
    }
}

void
input_matrix_d(double **A, size_t m, size_t n)
{
    int i, j;
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            scanf("%30lf", &A[j][i]);
        }
    }
}

void
free_matrix_d(double **A, size_t m, size_t n)
{
    int j;
    for (j = 0; j < n; ++j)
    {
        free(A[j]);
    }
    free(A);
}

void
matmul_dd(double **rop, double **A, double **B, size_t m, size_t n, size_t p)
{
    int i, j, k;
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < p; ++j)
        {
            rop[j][i] = 0;
            for (k = 0; k < n; ++k)
            {
                rop[j][i] += A[k][i] * B[j][k];
            }
        }
    }
}


void
gram_schmidt_d(double **rop, double **input, size_t m, size_t n)
{
    int i, j, k;

    if (rop != input)
    {
        copy_matrix_dd(rop, input, m, n);
    }
    for (k = 0; k < n; ++k)
    {
        double norm;
        norm = 0.0;
        for (i = 0; i < m; ++i)
        {
            norm += rop[k][i] * rop[k][i];
        }
        norm = sqrt(norm);
        if (norm > 1e-13)
        {
            for (i = 0; i < m; ++i)
            {
                rop[k][i] = rop[k][i] / norm;
            }
        }
        for (j = k + 1; j < n; ++j)
        {
            norm = col_product_dd(rop, rop, k, j, m);
            for (i = 0; i < m; ++i)
            {
                rop[j][i] = rop[j][i] - rop[k][i] * norm;
            }
        }
    }
}

void
qr_d(double **rop_q, double **rop_r, double **input, size_t m, size_t n)
{
    int i, j, k;

    if (rop_q != input)
    {
        copy_matrix_dd(rop_q, input, m, n);
    }
    for (k = 0; k < n; ++k)
    {
        double norm;
        norm = 0.0;
        for (i = 0; i < m; ++i)
        {
            norm += rop_q[k][i] * rop_q[k][i];
        }
        rop_r[k][k] = sqrt(norm);
        if (rop_r[k][k] > 1e-13)
        {
            for (i = 0; i < m; ++i)
            {
                rop_q[k][i] = rop_q[k][i] / rop_r[k][k];
            }
        }
        else
        {
            rop_r[k][k] = 0.0;
        }
        for (j = k + 1; j < n; ++j)
        {
            rop_r[j][k] = col_product_dd(rop_q, input, k, j, m);
            for (i = 0; i < m; ++i)
            {
                rop_q[j][i] = rop_q[j][i] - rop_q[k][i] * rop_r[j][k];
            }
        }
    }
}

int
main(void)
{
    int m, n;
    double **A, **B, **Q, **R, **QR;
    scanf("%20d", &m);
    scanf("%20d", &n);

    if (m < 1 || n < 1)
    {
        printf("The input matrix should be at least 1x1\n");
        return -1;
    }

    A = alloc_matrix_d(m, n);
    B = alloc_matrix_d(m, n);
    Q = alloc_matrix_d(m, n);
    R = alloc_matrix_d(n, n);
    QR = alloc_matrix_d(m, n);

    input_matrix_d(A, m, n);
    copy_matrix_dd(B, A, m, n);

#if PRINTOUTPUT == 1
    printf("Gram-Schmidt orthogonalization\n");
#endif
    gram_schmidt_d(B, A, m, n);
#if PRINTOUTPUT == 1
    print_matrix_d(B, m, n);
#endif

#if PRINTOUTPUT == 1
    printf("QR factorization\n");
#endif
    qr_d(Q, R, A, m, n);

#if PRINTOUTPUT == 1
    printf("Q = \n");
    print_matrix_d(Q, m, n);
    printf("R = \n");
    print_matrix_d(R, n, n);
#endif

    free_matrix_d(A, m, n);
    free_matrix_d(B, m, n);
    free_matrix_d(Q, m, n);
    free_matrix_d(R, n, n);
    free_matrix_d(QR, m, n);

    return 0;
}
