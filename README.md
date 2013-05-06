QR factorization
================

A program that calculates the Gram-Schmidt orthonormalization and the QR factorization (using the Gram-Schmidt orthogonalization, both of them as given by Golub).

The input should be in the following format, where `m` is the number of rows and `n` is the number of columns of the matrix:
```
m n
11 12 13 ... 1n
21 22 23 ... 2n
.
.
.
m1 m2 m3 ... mn
```

Some sample input files are available in `/test`.

Two versions are currently available: one "naive" using the algorithm directly, and one using a transposed internal representation to optimize cache access. The results are fairly interesting:

| Matrix size | Naive    | Optimized    |
| ----------- |---------:| ------------:|
| 100 x 100   | 0.013s   | 0.020s       |
| 500 x 500   | 2.067s   | 0.350s       |
| 1000 x 1000 | 39.306s  | 3.766s       |
