/*
The fifth-order, L-stable, fully-implicit Radau IIA method with 3 stages
*/

#ifndef ODE_RADAU_IIA_5_H_
#define ODE_RADAU_IIA_5_H_

#include "ode_base.h"
#include "ode_irk.h"
#include "ode_newton_bridge.h"

class OdeRadauIIA5;

class RadauIIA5Newton : public OdeNewtonIRK<OdeRadauIIA5> {
    public:
        RadauIIA5Newton (unsigned long neq, unsigned long nnew, OdeRadauIIA5 *integrator) : OdeNewtonIRK (neq, nnew, integrator) {};
    private:
        void f_Newton (double *x, double *y);
        void J_Newton (double *x, double **J);
};

class OdeRadauIIA5 : public OdeBase, private OdeIRK {
    //friends!
    friend class OdeNewtonBridge<OdeRadauIIA5>;
    friend class OdeNewtonIRK<OdeRadauIIA5>;

    public:
        OdeRadauIIA5 (unsigned long neq);
        ~OdeRadauIIA5 ();
        RadauIIA5Newton get_newton () { return(*newton_); }
    private:
        double **a;
        double *b;
        RadauIIA5Newton *newton_;
        void step_ (double dt);
};

#endif
