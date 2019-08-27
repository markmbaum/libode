//miscellaneous useful functions

#ifndef ODE_UTIL_H_
#define ODE_UTIL_H_

#include <cmath>

//!simple maximum of two doubles
double ode_max2 (double a, double b);

//!simple minimum of two doubles
double ode_min2 (double a, double b);

//!checks if two numbers are very close to each other
bool ode_is_close (double a, double b, double thresh);

#endif
