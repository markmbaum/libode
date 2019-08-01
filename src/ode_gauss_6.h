/*
The sixth-order, A-stable, fully-implicit Gauss-Legendre method with 3 stages
*/

#ifndef ODE_GAUSS_6_H_
#define ODE_GAUSS_6_H_

#include "ode_base.h"
#include "ode_irk.h"
#include "ode_newton_bridge.h"

class OdeGauss6;

class Gauss6Newton : public OdeNewtonIRK<OdeGauss6> {
    public:
        Gauss6Newton (unsigned long neq, unsigned long nnew, OdeGauss6 *integrator) : OdeNewtonIRK (neq, nnew, integrator) {};
    private:
        void f_Newton (double *x, double *y);
        void J_Newton (double *x, double **J);
};

class OdeGauss6 : public OdeBase, private OdeIRK {
    //friends!
    friend class OdeNewtonBridge<OdeGauss6>;
    friend class OdeNewtonIRK<OdeGauss6>;

    public:
        OdeGauss6 (unsigned long neq);
        ~OdeGauss6 ();
        Gauss6Newton get_newton () { return(*newton_); }
    private:
        double **a;
        double *b;
        Gauss6Newton *newton_;
        void step_ (double dt);
};

#endif
