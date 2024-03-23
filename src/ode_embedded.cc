//! \file ode_embedded.cc

#include "ode_embedded.h"

namespace ode {

OdeEmbedded::OdeEmbedded (unsigned long neq, bool need_jac, int lowerord) :
    OdeAdaptive (neq, need_jac) {

    //order of the LOWER order solution used for error estimation
    lowerord_ = lowerord;
    //variables for adaptive time stepping and step rejection
    facsafe_ = 0.9; //fraction to multiply facopt by to be cautious
    facmin_ = 1e-2; //minumum time step reduction factor
    facmax_ = 1e1; //maximum time step increase factor
    //lower order solution
    solemb_ = new double[neq];
}

//destructor
OdeEmbedded::~OdeEmbedded () {
    delete [] solemb_;
}

double OdeEmbedded::error (double abstol, double reltol) {

    //index
    unsigned long i;
    //storage for evaluating error terms
    double sc, d, g, err=0.0;

    //compute the Inf norm of the difference between higher and lower order steps
    //with a scaling factor that combines relative and absolute error
    for (i=0; i<neq_; i++) {
        d = fabs(solemb_[i] - sol_[i]);
        sc = abstol + reltol*fabs(sol_[i]);
        g = d/sc;
        if ( g > err ) err = g;
    }

    return(err);
}

double OdeEmbedded::facopt (double err) {

    //compute the "optimal" adjustment factor for the time step
    double fac = facsafe_*pow(1.0/err, 1.0/(lowerord_ + 1.0));
    //keep the factor within limits
    fac = ode_max2(facmin_, fac);
    fac = ode_min2(facmax_, fac);

    return(fac);
}

void OdeEmbedded::adapt (double abstol, double reltol) {

    //compute the error estimate
    double err = error(abstol, reltol);
    //determine if the step should be rejected
    isrej_ = ( err >= 1.0 ) ? true : false;
    //determine the next time step
    dtopt_ = facopt(err)*dt_;
}

bool OdeEmbedded::is_rejected () {
    //simply return the value set in the adapt() function
    return( isrej_ );
}

double OdeEmbedded::dt_adapt () {
    //simply return the value set in the adapt() function
    return( dtopt_ );
}

} // namespace ode 
