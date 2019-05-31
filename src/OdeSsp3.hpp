/*
This is a strong stability preserving method of order 3, from Shu and Osher:
    C. W. Shu and S. Osher, Effcient implementation of essentially nonoscillatory shock-capturing schemes, J. Comput. Phys., 77, 1988, pp. 439-471.
*/

#ifndef ODE_SSP3_HPP
#define ODE_SSP3_HPP

#include "OdeBaseRK.hpp"

class OdeSsp3 : public OdeBaseRK {
public:
    //constructor
    OdeSsp3 (unsigned long long neq_);
    //destructor
    ~OdeSsp3 ();
    //take a time step
    void step (double dt);
private:
    //vectors for time stepping
    double *k1, *k2, *k3;
    //coefficents of tableau
    double c2, a21,
           c3, a31, a32,
                b1,  b2, b3;
};


#endif
