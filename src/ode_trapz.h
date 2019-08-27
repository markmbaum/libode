
#ifndef ODE_TRAPZ_H_
#define ODE_TRAPZ_H_

#include "ode_base.h"
#include "ode_rk.h"
#include "ode_erk.h"

/*!
Second order, explicit trapezoidal rule
*/
class OdeTrapz : public OdeBase, private OdeRK, private OdeERK {

    public:
        //!constructs
        /*!
        \param[in] neq size of ODE sytem
        */
        OdeTrapz (unsigned long neq);

    private:
        //function for taking a single time step
        void step_ (double dt);
        //coefficents of tableau
        double c2, a21,
                    b1, b2;
};

#endif
