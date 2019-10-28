//! \file ode_irk.cc

#include "ode_irk.h"

OdeIRK::OdeIRK (unsigned long neq, int nk) {

    //store nk
    nk_ = nk;
    //make a single large vector for k values
    kall_ = new double[neq*nk];
    //assign the k_ values along the kall_ array
    k_ = new double*[nk];
    for (int i=0; i<nk; i++) k_[i] = kall_ + i*neq;
}

OdeIRK::~OdeIRK () {
    delete [] k_;
    delete [] kall_;
}
