//! \file ode_util.cc

#include "ode_util.h"

double ode_max2 (double a, double b) {
    if (a > b) return(a);
    return(b);
}

double ode_min2 (double a, double b) {
    if (a < b) return(a);
    return(b);
}

bool ode_is_close (double a, double b, double thresh) {

    //get magnitudes
    double absa = fabs(a);
    double absb = fabs(b);
    double absd = fabs(a - b);
    //check relative differnence against a threshold
    if ((absd/absa < thresh) && (absd/absb < thresh))
        return(true);
    //otherwise the numbers aren't close
    return(false);
}
