//! \file ode_erk.cc

#include "ode_erk.h"

namespace ode {

OdeERK::OdeERK (unsigned long neq) {

    //temporary solution vector
    soltemp_ = new double[neq];
}

//destructor
OdeERK::~OdeERK () {
    delete [] soltemp_;
}

} // namespace ode 
