#include "qr.h"

/*
    Modified Gram-Schmidt algorithm as given by Golub
*/

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
            norm += rop[i][k] * rop[i][k];
        }
        norm = sqrt(norm);
        if (norm > 1e-13)
        {
            for (i = 0; i < m; ++i)
            {
                rop[i][k] = rop[i][k] / norm;
            }
        }
        for (j = k + 1; j < n; ++j)
        {
            norm = col_product_dd(rop, rop, k, j, m);
            for (i = 0; i < m; ++i)
            {
                rop[i][j] = rop[i][j] - rop[i][k] * norm;
            }
        }
    }
}

void
gram_schmidt_md(mat_d_ptr rop, mat_d_ptr input)
{
    gram_schmidt_d(rop->matrix, input->matrix, input->rows, input->cols);
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
            norm += rop_q[i][k] * rop_q[i][k];
        }
        rop_r[k][k] = sqrt(norm);
        if (rop_r[k][k] > 1e-13)
        {
            for (i = 0; i < m; ++i)
            {
                rop_q[i][k] = rop_q[i][k] / rop_r[k][k];
            }
        }
        else
        {
            rop_r[k][k] = 0.0;
        }
        for (j = k + 1; j < n; ++j)
        {
            rop_r[k][j] = col_product_dd(rop_q, input, k, j, m);
            for (i = 0; i < m; ++i)
            {
                rop_q[i][j] = rop_q[i][j] - rop_q[i][k] * rop_r[k][j];
            }
        }
    }
}

void
qr_md(mat_d_ptr rop_q, mat_d_ptr rop_r, mat_d_ptr input)
{
    qr_d(rop_q->matrix, rop_r->matrix, input->matrix, input->rows,
         input->cols);
}

int
main(void)
{
    int m, n;
    double **A, **B, **Q, **R, **QR;
    mat_d_t A2, B2, Q2, R2;
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

    mat_d_init(A2, m, n);
    mat_d_init(B2, m, n);
    mat_d_init(Q2, m, n);
    mat_d_init(R2, n, n);

    input_matrix_d(A, m, n);
    copy_matrix_dd(B, A, m, n);

    mat_d_set_d(A2, A);
    mat_d_set_md(B2, A2);

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

#if PRINTOUTPUT == 1
    printf("New implementation with structs\n");
    printf("Gram-Schmidt orthogonalization\n");
    gram_schmidt_md(B2, A2);
    mat_d_print(B2);

    printf("QR factorization\n");
    qr_md(Q2, R2, A2);
    printf("Q = \n");
    mat_d_print(Q2);
    printf("R = \n");
    mat_d_print(R2);
#endif

    mat_d_clear(A2);
    mat_d_clear(B2);
    mat_d_clear(Q2);
    mat_d_clear(R2);

    return 0;
}
