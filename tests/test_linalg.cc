#include <cmath>

#include "ode_linalg.h"

using namespace ode;

int main () {

    //size of system
    int n = 4;
    //matrices
    double **A = new double*[n];
    for (int i=0; i<n; i++) A[i] = new double[n];
    double **T = new double*[3];
    for (int i=0; i<3; i++) T[i] = new double[n];
    //permutation array
    int *p = new int[n];
    //right side
    double *b = new double[n];
    //temp array
    double *temp = new double[n];
    //solution
    double *x = new double[n];

    //set matrix values
    A[0][0] = 1.0;  A[0][1] = 2.0;  A[0][2] = 3.0; A[0][3] = 10;
    A[1][0] = 5.0;  A[1][1] = -2.;  A[1][2] = 6.0; A[1][3] = -10;
    A[2][0] = 7.0;  A[2][1] = 8.0;  A[2][2] = -1.; A[2][3] = 0.01;
    A[3][0] = 10.0;  A[3][1] = -4.;  A[3][2] = 12.; A[3][3] = -19;

    //fill tridiagonal matrix (before A is LU decomposed in-place)
    T[0][0] = NAN;
    for (int i=0; i<n-1; i++) T[0][i+1] = A[i][i+1];
    for (int i=0; i<n;   i++) T[1][i]   = A[i][i];
    for (int i=0; i<n-1; i++) T[2][i]   = A[i+1][i];
    T[2][n-1] = NAN;

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
    //should be [-22.47104762,  22.09780952,  16.42514286,  -6]
    printf("x =\n");
    for (int i=0; i<n; i++) {
        printf("%16g\n", x[i]);
    }

    //print tridiagonal version of the system
    printf("T =\n");
    for (int i=0; i<3; i++) {
        for (int j=0; j<n-1; j++) {
            printf("%8g ", T[i][j]);
        }
        printf("%8g\n", T[i][n-1]);
    }

    //solve
    ode_solve_tridiag(T, b, temp, n, x);

    //print the solution
    //should be [12.92215219,  -0.9610761 , -10.75548553,  -6.68767507]
    printf("x =\n");
    for (int i=0; i<n; i++) {
        printf("%16g\n", x[i]);
    }

    for (int i=0; i<n; i++) delete [] A[i];
    delete [] A;
    for (int i=0; i<3; i++) delete [] T[i];
    delete [] T;
    delete [] p;
    delete [] b;
    delete [] x;
    delete [] temp;
}
