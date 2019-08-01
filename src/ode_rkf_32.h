/*
2nd and 3rd order solver developed by Fehlberg
*/

#ifndef ODE_RKF32_H_
#define ODE_RKF32_H_

#include "ode_embedded.h"
#include "ode_rk.h"
#include "ode_erk.h"

class OdeRKF32 : public OdeEmbedded, private OdeRK, private OdeERK {

    public:
        //constructor
        OdeRKF32 (unsigned long neq);

    private:
        //function for taking a single time step
        void step_ (double dt);
        //coefficents of tableau
        double c2, a21,
               c3, a31, a32,
                    b1,  b2, b3,
                    d1,  d2, d3;
};

#endif
