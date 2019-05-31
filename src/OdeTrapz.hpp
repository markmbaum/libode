/*
OdeTrapz is an adaptive version of the trapezoidal method, using the first order Euler solution for error estimation. The solver can be used non-adaptively.
*/

#ifndef ODE_TRAPZ_HPP
#define ODE_TRAPZ_HPP

#include "OdeBaseARK.hpp"

class OdeTrapz : public OdeBaseARK {
public:
    //constructor
    OdeTrapz (unsigned long long neq_);
    //destructor
    ~OdeTrapz ();
    //function for taking a single time step
    void step (double dt);
private:
    //vectors for stepping
    double *k1, *k2;
    //coefficents of tableau
    double c2, a21,
                b1, b2;

};

#endif
