/*
These classes provides some variables of convenience when connecting the
Newton base (which is left as a generic Newton solver for potential reuse
elsewhere) and the systems required for implicit methods.
*/

#ifndef ODE_NEWTON_BRIDGE_H_
#define ODE_NEWTON_BRIDGE_H_

#include "ode_newton.h"

//general connector
template <class T>
class OdeNewtonBridge : public OdeNewton {

    public:
        //constructor
        OdeNewtonBridge (unsigned long neq, unsigned long nnew, T *integrator)
            : OdeNewton (nnew) {

            //make the Newton iteration limit quite high
            set_iter_Newton(10000);

            //store pointers to the ODE solver object and its members
            integrator_ = integrator;
            sol_ = integrator_->sol_;
            Jac_ = integrator_->Jac_;
            dt_ = &(integrator_->dt_);
            //store ODE system size
            neq_ = integrator_->neq_;
            //store Newton system size;
            nnew_ = nnew;

            //temporary values for evaluation of ODE fun within Newton fun
            ftemp_ = new double[neq];
            //temporary solution values
            soltemp_ = new double[neq];
        }
        //destructor
        ~OdeNewtonBridge () {
            delete [] ftemp_;
            delete [] soltemp_;
        }

    protected:
        //system sizes
        unsigned long neq_, nnew_;
        //storage of a pointer to the solver class
        T *integrator_;
        //pointer to the solver's solution vector
        double *sol_;
        //pointer to the solver's Jacobian matrix
        double **Jac_;
        //pointer to time step member
        double *dt_;

        //temporary values for evaluation of Newton function
        double *ftemp_;
        //temporary solution values
        double *soltemp_;

        //wrappers around solver functions
        void fun(double *solin, double *fout) { integrator_->ode_fun_(solin, fout); }
        void jac(double *solin, double **Jout) { integrator_->ode_jac_(solin, Jout); }
};

/*
This class extends the Bridge class by storing pointers to the integrator class's
k vectors and tableau/coefficients, which are needed for newton iterations
*/
template <class T>
class OdeNewtonIRK : public OdeNewtonBridge<T> {

    public:
        OdeNewtonIRK (unsigned long neq, int nk, T *integrator) :
            OdeNewtonBridge<T> (neq, neq*nk, integrator) {

            //number of stages or k vectors
            nk_ = nk;
            //pointer to stage k values
            k_ = integrator->k_;
            //pointers to tableau coefficients
            a = integrator->a;
            b = integrator->b;
        }

    protected:
        //number of stages or k vectors
        int nk_;
        //pointer to the stage slopes of RK methods
        double **k_;
        //pointers to tableau coefficients
        double **a;
        double *b;
};

/*
This class extends the Bridge class by storing pointers to the integrator class's
k vectors and tableau/coefficients, which are needed for newton iterations. It
also stores an ik_ integer indicating which stage of the SDIRK scheme is being
solved for and the gamma value of the scheme (the diagonal of the a matrix).
*/
template <class T>
class OdeNewtonSDIRK : public OdeNewtonBridge<T> {

    public:
        OdeNewtonSDIRK (unsigned long neq, T *integrator) :
            OdeNewtonBridge<T> (neq, neq, integrator) {

            //pointer to stage k values
            k_ = integrator->k_;
            //pointers to tableau coefficients
            gam = integrator->gam;
            a = integrator->a;
            b = integrator->b;
            //index of k vector being solved for
            ik_ = 0;
        }

        void set_ik (int ik) { ik_ = ik; }

    protected:
        //pointer to the stage slopes of RK methods
        double **k_;
        //pointers to tableau coefficients
        double gam;
        double **a;
        double *b;
        //the index of the k vector being solved for
        int ik_;
};


#endif
