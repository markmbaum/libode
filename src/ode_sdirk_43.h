/*
L-stable 4/3 SDIRK pair from section IV.6 of:
    Wanner, Gerhard, and Ernst Hairer. Solving ordinary differential equations II. Springer Berlin Heidelberg, 1996.
*/

#ifndef ODE_SDIRK_43_H_
#define ODE_SDIRK_43_H_

#include "ode_embedded.h"
#include "ode_rk.h"
#include "ode_newton_bridge.h"

class OdeSDIRK43;

class SDIRK43Newton : public OdeNewtonSDIRK<OdeSDIRK43> {
    public:
        SDIRK43Newton (unsigned long neq, OdeSDIRK43 *integrator) : OdeNewtonSDIRK (neq, integrator) {};
    private:
        void f_Newton (double *x, double *y);
        void J_Newton (double *x, double **J);
};

class OdeSDIRK43 : public OdeEmbedded, private OdeRK {
    //friends!
    friend class OdeNewtonBridge<OdeSDIRK43>;
    friend class OdeNewtonSDIRK<OdeSDIRK43>;

    public:
        OdeSDIRK43 (unsigned long neq);
        ~OdeSDIRK43 ();
        SDIRK43Newton get_newton () { return(*newton_); }
    private:
        double gam;
        double **a;
        double *b;
        double *d;
        SDIRK43Newton *newton_;
        void step_ (double dt);
};

#endif
