/*
OdeMidpnt is an adaptive version of the midpoint method, using the first order Euler solution for error estimation and advancing with the second order midpoint solution. It can be used without adaptive time stepping.
*/

#ifndef ODE_MIDPNT_HPP
#define ODE_MIDPNT_HPP

#include "OdeBaseARK.hpp"

class OdeMidpnt : public OdeBaseARK {
public:
    //constructor
    OdeMidpnt (unsigned long long neq_);
    //destructor
    ~OdeMidpnt ();
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
