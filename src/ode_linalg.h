#ifndef ODE_LINALG_H_
#define ODE_LINALG_H_

#include <cmath>
#include <cstdio>
#include <cstdlib>

/*performs forward substitution to solve a lower triangular matrix equation
where the elements on the main diagonal of the lower triangular matrix are \
assumed to all be 1
input:
    L - n by n lower triangular matrix
    b - right hand side of matrix equation U*out = b
    p - array of permutation indices from crout_LU
    n - dimensions of A and length of b
output:
    out - array of length n, result of substitution*/
void ode_crout_forw_sub (double **L, double *b, int *p, int n, double *out);

/*performs backward substitution to solve an upper triangular matrix equation
input:
    U - n by n upper triangular matrix
    b - right hand side of matrix equation U*out = b
    n - dimensions of A and length of b
output:
    out - array of length n, result of substitution*/
void ode_back_sub (double **U, double *b, int n, double *out);

/*performs LU decomposition with partial pivoting on the matrix A and sets
permutation indices in the array p, following the outline in the book:
    Demmel, James W. Applied numerical linear algebra. Vol. 56. Siam, 1997.
input:
    A - matrix of size n to decompose (overwritten)
    n - size of square matrix
output:
    A - the factorized matrix overwrites the input matrix A
    p - integer array loaded with permutation indices*/
void ode_crout_LU (double **A, int n, int *p);

/*solves a matrix equation where the matrix has already be crout LU decomposed
input:
    LU - crout LU decomposed matrix
    p - permutation indices
    b - right-hand side of matrix equation
    n - size of matrix and arrays
output:
    out - solution to linear system*/
void ode_solve_LU (double **LU, int *p, double *b, int n, double *out);

/*solves a matrix equation A x = b with crout LU decomposition
input:
    A - matrix, is overwritten

    b - right-hand side of matrix equation
    n - size of matrix and arrays
output:
    out - solution to linear system*/
void ode_solve_A (double **A, double *b, int n, double *out);

#endif
