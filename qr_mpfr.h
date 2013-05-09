#ifndef QR_MPFR_H
#define QR_MPFR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>

#define ROUNDINGMODE 0

/*
	MPFR-based matrix struct
*/
typedef struct
{
    size_t rows;
    size_t cols;
    mpfr_prec_t prec;
    mpfr_t *matrix;
} __struct_mat_fr;

typedef __struct_mat_fr mat_fr_t[1];
typedef __struct_mat_fr *mat_fr_ptr;

/*
	Main functions
*/
void gram_schmidt_mat_fr(mat_fr_ptr rop, mat_fr_ptr input, mpfr_rnd_t rnd);
/*
	Utility functions
*/
mpfr_t * alloc_matrix_fr(size_t m, size_t n, mpfr_prec_t prec);
void mat_fr_init(mat_fr_ptr A, size_t m, size_t n, mpfr_prec_t prec);
int mat_fr_set_fr(mat_fr_ptr A, mpfr_ptr op, size_t row, size_t col, mpfr_rnd_t rnd);
int mat_fr_set_si(mat_fr_ptr A, long int op, size_t row, size_t col, mpfr_rnd_t rnd);
int mat_fr_set_d(mat_fr_ptr A, double op, size_t row, size_t col, mpfr_rnd_t rnd);
mpfr_ptr mat_fr_get(mat_fr_ptr A, size_t row, size_t col);
void mat_fr_set(mat_fr_ptr rop, mat_fr_ptr A, mpfr_rnd_t rnd);
void mat_fr_clear(mat_fr_ptr A);
void mat_fr_col_prod(mpfr_ptr rop, mat_fr_ptr A, mat_fr_ptr B, size_t a, size_t b, size_t m, mpfr_rnd_t rnd);
void mat_fr_input_d(mat_fr_ptr A, mpfr_rnd_t rnd);
void mat_fr_print(mat_fr_ptr A);

#endif
