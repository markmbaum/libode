#include "ode_misc.hpp"

double ode_max2 (double a, double b) {
    if (a > b) return(a);
    return(b);
}

double ode_min2 (double a, double b) {
    if (a < b) return(a);
    return(b);
}

void ode_write_double (char const *fn, double *a, long size) {

    ode_check_file_write(fn);
    FILE *ofile;
    ofile = fopen(fn, "wb");
    fwrite(a, sizeof(double), size, ofile);
    fclose(ofile);
}


void ode_check_file_write (const char *fn) {
    FILE* ofile;
    ofile = fopen(fn, "wb");
    if (ofile == NULL) {
        fputs ("cannot open file\n", stderr);
        exit(EXIT_FAILURE);
    }
    fclose(ofile);
}

bool ode_is_close (double a, double b, double thresh) {
    //check equality
    if (a == b) return(true);
    //check signs
    if ((a > 0 && b > 0) || (a < 0 && b < 0)) {
        //remove signs
        double absa = fabs(a);
        double absb = fabs(b);
        double absd = fabs(a - b);
        //check relative differnence against a threshold
        if ((absd/absa < thresh) && (absd/absb < thresh)) return(true);
    }
    //otherwise the numbers aren't close
    return(false);
}
