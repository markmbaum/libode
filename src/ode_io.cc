#include "ode_io.h"

void ode_check_write (const char *fn) {

    //try to open the file
    FILE *ofile = fopen(fn, "wb");
    //raise an error if needed
    if (ofile == NULL) ode_print_exit("cannot open file, ensure the requisite directories exist (like the output directory)");
    //close
    fclose(ofile);
}

void ode_print_exit (const char *msg) {

    //print the given message
    printf("\nODE FAILURE: %s\n\n", msg);
    //cancel the program
    exit(EXIT_FAILURE);
}

std::string int_to_string (long i) {

    // https://www.geeksforgeeks.org/converting-string-to-number-and-vice-versa-in-c/
    std::ostringstream str;
    str << i;
    std::string out = str.str();
    return(out);
}
