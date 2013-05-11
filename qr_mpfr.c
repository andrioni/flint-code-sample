#include "qr_mpfr.h"

mpfr_t *
alloc_matrix_fr(size_t m, size_t n, mpfr_prec_t prec)
{
    int i;
    mpfr_t *A;
    A = malloc(m*n*sizeof(mpfr_t));
    for (i = 0; i < m*n; ++i)
    {
        mpfr_init2(A[i], prec);
    }
    return A;
}

void
mat_fr_init(mat_fr_ptr A, size_t m, size_t n, mpfr_prec_t prec)
{
    A->rows = m;
    A->cols = n;
    A->prec = prec;
    A->matrix = alloc_matrix_fr(m, n, prec);
}

int
mat_fr_set_fr(mat_fr_ptr A, mpfr_ptr op, size_t row, size_t col, mpfr_rnd_t rnd)
{
    return mpfr_set(A->matrix[row + col * A->rows], op, rnd);
}

int
mat_fr_set_si(mat_fr_ptr A, long int op, size_t row, size_t col, mpfr_rnd_t rnd)
{
    return mpfr_set_si(A->matrix[row + col * A->rows], op, rnd);
}

int
mat_fr_set_d(mat_fr_ptr A, double op, size_t row, size_t col, mpfr_rnd_t rnd)
{
    return mpfr_set_d(A->matrix[row + col * A->rows], op, rnd);
}

mpfr_ptr
mat_fr_get(mat_fr_ptr A, size_t row, size_t col)
{
    return A->matrix[row + col * A->rows];
}

void
mat_fr_set(mat_fr_ptr rop, mat_fr_ptr A, mpfr_rnd_t rnd)
{
    int i, j;
    for (i = 0; i < rop->rows; ++i)
    {
        for (j = 0; j < rop->cols; ++j)
        {
            mpfr_set(mat_fr_get(rop, i, j), mat_fr_get(A, i, j), rnd);
        }
    }
}

void
mat_fr_clear(mat_fr_ptr A)
{
    int i;
    for (i = 0; i < A->rows * A->cols; ++i)
    {
        mpfr_clear(A->matrix[i]);
    }
}

void
mat_fr_col_prod(mpfr_ptr rop, mat_fr_ptr A, mat_fr_ptr B, size_t a, size_t b, size_t m, mpfr_rnd_t rnd)
{
    int i;
    mpfr_t temp;
    mpfr_init2(temp, mpfr_get_prec(rop));
    mpfr_set_si(temp, 0, rnd);
    mpfr_set_si(rop, 0, rnd);
    for (i = 0; i < m; ++i)
    {
        mpfr_mul(temp, mat_fr_get(A, i, a), mat_fr_get(B, i, b), rnd);
        mpfr_add(rop, rop, temp, rnd);
    }
    mpfr_clear(temp);
}

void
mat_fr_input_d(mat_fr_ptr A, mpfr_rnd_t rnd)
{
    int i, j;
    double aux;
    for (i = 0; i < A->rows; ++i)
    {
        for (j = 0; j < A->cols; ++j)
        {
            scanf("%30lf", &aux);
            mat_fr_set_d(A, aux, i, j, rnd);
        }
    }
}

void
mat_fr_print(mat_fr_ptr A)
{
    int i, j;
    for (i = 0; i < A->rows; ++i)
    {
        for (j = 0; j < A->cols; ++j)
        {
            mpfr_printf("%20.10Rf ", mat_fr_get(A, i, j));
        }
        printf("\n");
    }
}

void
mpfrdp(mpfr_ptr op)
{
    mpfr_printf("%30.15Rf\n", op);
}


void
gram_schmidt_mat_fr(mat_fr_ptr rop, mat_fr_ptr input, mpfr_rnd_t rnd)
{
    int i, j, k, m, n;
    mpfr_t temp, norm;

    m = rop->rows;
    n = rop->cols;
    mpfr_init2(temp, rop->prec);
    mpfr_init2(norm, rop->prec);

    if (rop != input)
    {
        mat_fr_set(rop, input, rnd);
    }
    for (k = 0; k < n; ++k)
    {
        mpfr_set_si(norm, 0, rnd);
        for (i = 0; i < m; ++i)
        {
            mpfr_mul(temp, mat_fr_get(rop, i, k), mat_fr_get(rop, i, k), rnd);
            mpfr_add(norm, norm, temp, rnd);
        }
        mpfr_sqrt(norm, norm, rnd);
        mpfr_set_si(temp, 0, rnd);
        if (!mpfr_equal_p(norm, temp))
        {
            for (i = 0; i < m; ++i)
            {
                mpfr_div(mat_fr_get(rop, i, k), mat_fr_get(rop, i, k), norm, rnd);
            }
        }
        for (j = k + 1; j < n; ++j)
        {
            mat_fr_col_prod(norm, rop, rop, k, j, m, rnd);
            for (i = 0; i < m; ++i)
            {
                mpfr_mul(temp, mat_fr_get(rop, i, k), norm, rnd);
                mpfr_sub(mat_fr_get(rop, i, j), mat_fr_get(rop, i, j), temp, rnd);
            }
        }
    }

    mpfr_clear(temp);
    mpfr_clear(norm);
}

int
main(void)
{
    int m, n;
    mat_fr_t A;

    scanf("%20d", &m);
    scanf("%20d", &n);

    mat_fr_init(A, m, n, 256);
    mat_fr_input_d(A, ROUNDINGMODE);
    gram_schmidt_mat_fr(A, A, ROUNDINGMODE);

    printf("Gram-Schmidt orthogonalization\n");
    mat_fr_print(A);
    mat_fr_clear(A);
    return 0;
}
