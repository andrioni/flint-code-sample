#include "qr.h"

double **
alloc_matrix_d(size_t m, size_t n)
{
    int i;
    double **A;
    A = calloc(m, sizeof(double *));
    for (i = 0; i < m; ++i)
    {
        A[i] = calloc(n, sizeof(double));
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
        result += A[i][a] * B[i][b];
    }
    return result;
}

void
copy_matrix_dd(double **rop, double **A, size_t m, size_t n)
{
    int i, j;
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            rop[i][j] = A[i][j];
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
            printf("%10.5f ", matrix[i][j]);
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
            scanf("%30lf", &A[i][j]);
        }
    }
}

void
free_matrix_d(double **A, size_t m, size_t n)
{
    int i;
    for (i = 0; i < m; ++i)
    {
        free(A[i]);
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
            rop[i][j] = 0;
            for (k = 0; k < n; ++k)
            {
                rop[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void
mat_d_init(mat_d_ptr A, size_t m, size_t n)
{
    A->rows = m;
    A->cols = n;
    A->matrix = alloc_matrix_d(m, n);
}

void
mat_d_input(mat_d_ptr A)
{
    input_matrix_d(A->matrix, A->rows, A->cols);
}

void
mat_d_clear(mat_d_ptr A)
{
    free_matrix_d(A->matrix, A->rows, A->cols);
}

void
mat_d_print(mat_d_ptr A)
{
    print_matrix_d(A->matrix, A->rows, A->cols);
}

void
mat_d_set_d(mat_d_ptr rop, double **inp)
{
    copy_matrix_dd(rop->matrix, inp, rop->rows, rop->cols);
}

void
mat_d_set_md(mat_d_ptr rop, mat_d_ptr inp)
{
    mat_d_clear(rop);
    mat_d_init(rop, inp->rows, inp->cols);
    copy_matrix_dd(rop->matrix, inp->matrix, inp->rows, inp->cols);
}
