#ifndef QR_H
#define QR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef PRINTOUTPUT
#define PRINTOUTPUT 1
#endif

/*
	Main functions
*/
void gram_schmidt_d(double **rop, double **input, size_t m, size_t n);
void qr_d(double **rop_q, double **rop_r, double **input, size_t m, size_t n);

/*
	Utility functions
*/
double **alloc_matrix_d(size_t m, size_t n);
void copy_matrix_dd(double **rop, double **A, size_t m, size_t n);
double col_product_dd(double **A, double **B, size_t a, size_t b, size_t m);
void print_matrix_d(double **matrix, size_t m, size_t n);
void input_matrix_d(double **A, size_t m, size_t n);
void free_matrix_d(double **A, size_t m);
void matmul_dd(double **rop, double **A, double **B, size_t m, size_t n,
               size_t p);

#endif
