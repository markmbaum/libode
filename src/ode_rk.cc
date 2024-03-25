//! \file ode_rk.cc

#include "ode_rk.h"

namespace ode {

OdeRK::OdeRK (unsigned long neq, int nk) {

    //store number of stages
    nk_ = nk;
    //storage for k values of stages
    k_ = new double*[nk];
    for (int i=0; i<nk; i++) k_[i] = new double[neq];
}

//destructor
OdeRK::~OdeRK () {
    for (int i=0; i<nk_; i++) delete [] k_[i];
    delete [] k_;
}

} // namespace ode
