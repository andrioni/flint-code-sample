#ifndef QR_H
#define QR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>

#ifndef PRINTOUTPUT
#define PRINTOUTPUT 1
#endif

/*
	Double-based matrix struct
*/
typedef struct
{
    size_t rows;
    size_t cols;
    double **matrix;
} __struct_mat_d;

typedef __struct_mat_d mat_d_t[1];
typedef __struct_mat_d *mat_d_ptr;

/*
	Main functions
*/
void gram_schmidt_d(double **rop, double **input, size_t m, size_t n);
void qr_d(double **rop_q, double **rop_r, double **input, size_t m, size_t n);
void gram_schmidt_d(double **rop, double **input, size_t m, size_t n);
void gram_schmidt_md(mat_d_ptr rop, mat_d_ptr input);

/*
	Utility functions
*/
double **alloc_matrix_d(size_t m, size_t n);
void copy_matrix_dd(double **rop, double **A, size_t m, size_t n);
double col_product_dd(double **A, double **B, size_t a, size_t b, size_t m);
void print_matrix_d(double **matrix, size_t m, size_t n);
void input_matrix_d(double **A, size_t m, size_t n);
void free_matrix_d(double **A, size_t m, size_t n);
void matmul_dd(double **rop, double **A, double **B, size_t m, size_t n,
               size_t p);
void mat_d_init(mat_d_ptr A, size_t m, size_t n);
void mat_d_input(mat_d_ptr A);
void mat_d_clear(mat_d_ptr A);
void mat_d_print(mat_d_ptr A);
void mat_d_set_d(mat_d_ptr rop, double **inp);
void mat_d_set_md(mat_d_ptr rop, mat_d_ptr inp);

#endif
