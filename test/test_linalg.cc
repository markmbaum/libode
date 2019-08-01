#include "ode_linalg.h"

int main () {

    //size of system
    int n = 4;
    //matrix
    double **A = new double*[n];
    for (int i=0; i<n; i++) A[i] = new double[n];
    //permutation array
    int *p = new int[n];
    //right side
    double *b = new double[n];
    //solution
    double *x = new double[n];

    //set matrix values
    A[0][0] = 1.0;  A[0][1] = 2.0;  A[0][2] = 3.0; A[0][3] = 10;
    A[1][0] = 5.0;  A[1][1] = -2.;  A[1][2] = 6.0; A[1][3] = -10;
    A[2][0] = 7.0;  A[2][1] = 8.0;  A[2][2] = -1.; A[2][3] = 0.01;
    A[3][0] = 10.0;  A[3][1] = -4.;  A[3][2] = 12.; A[3][3] = -19;

    //set right hand side
    b[0] = 11.;
    b[1] = 2.0;
    b[2] = 3.0;
    b[3] = -2.;

    //print the input elements
    printf("A =\n");
    for (int i=0; i<n; i++) {
        for (int j=0; j<n-1; j++) {
            printf("%8g ", A[i][j]);
        }
        printf("%8g\n", A[i][n-1]);
    }
    printf("b =\n");
    for (int i=0; i<n; i++) {
        printf("%8g\n", b[i]);
    }

    //decompose the matrix in-place
    ode_crout_LU(A, n, p);

    //solve the system
    ode_solve_LU(A, p, b, n, x);

    //print the solution
    printf("x =\n");
    for (int i=0; i<n; i++) {
        printf("%16g\n", x[i]);
    }

    for (int i=0; i<n; i++) delete [] A[i];
    delete [] A;
    delete [] p;
    delete [] b;
    delete [] x;
}
