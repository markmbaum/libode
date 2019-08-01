//input and output functions

#ifndef ODE_IO_H_
#define ODE_IO_H_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>

/*check if a file can be written (opened in write mode)
input:
    fn - desired file name*/
void ode_check_write (const char *fn);

/*write an array to a binary file
inputs:
    fn - file name or path
    a - array of values to write
    n - size of array or number of elements to write*/
template<class T>
void ode_write (char const *fn, T *a, unsigned long n) {
    //check that the file can be written
    ode_check_write(fn);
    //open the file
    FILE *ofile = fopen(fn, "wb");
    //write values
    fwrite(a, sizeof(T), n, ofile);
    //close the file
    fclose(ofile);
}

/*print a message and exit with failure
input:
    msg - the message to print before exiting the program*/
void ode_print_exit (const char *msg);

/*convert an integer into a string*/
std::string int_to_string (long i);

#endif
