//miscellaneous useful functions

#ifndef MISC_HPP
#define MISC_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>

double ode_max2 (double a, double b);
double ode_min2 (double a, double b);
void ode_write_double (char const *fn, double *a, long size);
void ode_check_file_write (const char *fn);
bool ode_is_close (double a, double b, double thresh);

#endif
