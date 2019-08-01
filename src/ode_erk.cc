#include "ode_erk.h"

OdeERK::OdeERK (unsigned long neq) {

    //temporary solution vector
    soltemp_ = new double[neq];
}

//destructor
OdeERK::~OdeERK () {
    delete [] soltemp_;
}
