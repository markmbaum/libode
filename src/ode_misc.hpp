//miscellaneous useful functions

#ifndef MISC_HPP
#define MISC_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>

//simple maximum of two doubles
double ode_max2 (double a, double b);

//simple minimum of two doubles
double ode_min2 (double a, double b);

//write an array to a binary file
void ode_write_double (char const *fn, double *a, long size);

//check of a file can be written
void ode_check_file_write (const char *fn);

//check of two numbers are very close to each other
bool ode_is_close (double a, double b, double thresh);

#endif
