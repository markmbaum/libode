/*
The sixth-order, L-stable, fully-implicit Lobatto IIIC method with 4 stages
*/

#ifndef ODE_LOBATTO_IIIC_6_H_
#define ODE_LOBATTO_IIIC_6_H_

#include "ode_base.h"
#include "ode_irk.h"
#include "ode_newton_bridge.h"

class OdeLobattoIIIC6;

class LobattoIIIC6Newton : public OdeNewtonIRK<OdeLobattoIIIC6> {
    public:
        LobattoIIIC6Newton (unsigned long neq, unsigned long nnew, OdeLobattoIIIC6 *integrator) : OdeNewtonIRK (neq, nnew, integrator) {};
    private:
        void f_Newton (double *x, double *y);
        void J_Newton (double *x, double **J);
};

class OdeLobattoIIIC6 : public OdeBase, private OdeIRK {
    //friends!
    friend class OdeNewtonBridge<OdeLobattoIIIC6>;
    friend class OdeNewtonIRK<OdeLobattoIIIC6>;

    public:
        OdeLobattoIIIC6 (unsigned long neq);
        ~OdeLobattoIIIC6 ();
        LobattoIIIC6Newton get_newton () { return(*newton_); }
    private:
        double **a;
        double *b;
        LobattoIIIC6Newton *newton_;
        void step_ (double dt);
};

#endif
