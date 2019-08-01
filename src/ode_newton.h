/*
This class implements Newton's method for nonlinear systems of equations.
The virtual functions F_Newton and Jac_Newton allow a derived class to implement
the system of equations and its Jacobian matrix.
*/

#ifndef ODE_NEWTON_H_
#define ODE_NEWTON_H_

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "ode_linalg.h"

class OdeNewton {

    public:
        //constructor
        OdeNewton (unsigned long n);
        //destructor
        ~OdeNewton ();

        //getters and setters
        unsigned long get_n () { return(n_); }
        double get_tol_Newton () { return(tol_Newton_); }
        unsigned long get_iter_Newton () { return(iter_Newton_); }
        int  get_iJLU () { return(iJLU_); }
        unsigned long get_nJLU () { return(nJLU_); }
        unsigned long get_n_solve_LU () { return(n_solve_LU_); }
        bool get_modified () { return(modified_); }
        bool get_ignore_JLU () { return(ignore_JLU_); }

        void set_tol_Newton (double tol_Newton) { tol_Newton_ = tol_Newton; }
        void set_iter_Newton (unsigned long iter_Newton) { iter_Newton_ = iter_Newton; }
        void set_iJLU (int iJLU) { iJLU_ = iJLU; }
        void set_modified (bool modified) { modified_ = modified; }
        void set_ignore_JLU (bool ignore_JLU) { ignore_JLU_ = ignore_JLU; }

        /*solve for a root of the function f_Newton using the jacobian J_Newton,
        both of which must be implemented in derived classes
        input:
            x - initial guess for the root
        output:
            x - final values for root
        returns:
            suc - success code, which is zero for success and
                    1 - too many iterations
                    2 - solution contains nan(s)
                    3 - solution contains inf(s)*/
        int solve_Newton (double *x);

    protected:
        //evaluates the function being zeroed
        virtual void f_Newton (double *x, double *f) = 0;
        //evaluates the Jacobian matrix of the function being zeroed
        virtual void J_Newton (double *x, double **J) = 0;

    private:
        //size of system
        unsigned long n_;
        //whether to use only a single evaluation and LU decomposition of Jac
        bool modified_;
        //whether to do no JLU updates at all
        bool ignore_JLU_;
        //error tolerance
        double tol_Newton_;
        //iteration tolerance
        unsigned long iter_Newton_;
        //number of iterations after which the Jacobian is updated and redecomposed
        int iJLU_;
        //number of times Jac is updated
        unsigned long nJLU_;
        //number of times the LU decomposed matrix is used to solve a matrix eq
        unsigned long n_solve_LU_;
        //iteration interval for checking solution integrity
        unsigned long icheck_;
        //current evaluation of F
        double *f_;
        //current evaluation of J
        double **J_;
        //update array
        double *delx_;
        //permutation array
        int *p_;
        //finds the infinity norm of the update vector delx and of y
        void err (double *errx, double *erry);
        //recomputes the Jacobian and crout LU decomposes it
        void JLU (double *x);
        //wrapper around LU solving routine for counting solves
        void solve_LU_(double **LU, int *p, double *b, int n, double *out);
        //function for checking solution integrity (nans and infs)
        int check_integrity (double *x);
};

#endif
