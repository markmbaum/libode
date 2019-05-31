#include "OdeBase.hpp"

OdeBase::OdeBase (unsigned long long neq_) {

    //number of equations/variables in ode system
    neq = neq_;
    //time starts at zero unless modified
    t = 0.0;
    //time step
    dt = NAN;
    //number of time steps, must be incremented in step()
    nstep = 0;
    //number of function evaluations, must be incremented in step()
    neval = 0;
    //interval of solution integrity checks
    ncheck = 100;
    //array of solution values for a given time
    sol = new double[neq];
}

OdeBase::~OdeBase () {
    delete [] sol;
}

void OdeBase::snap (std::string dirout, long *isnap, double *tsnap,
    double tend, int snaps, std::vector<double> *tout) {

    //output filename
    std::string fnout = dirout + "/snap_" + std::to_string(*isnap);
    //write output
    ode_write_double(fnout.data(), sol, neq);
    //progress report
    printf("    snap %li written\n", *isnap);
    //do any extra snapping
    snap_extra(dirout, *isnap, *tsnap);

    //update
    tout->push_back(t);
    *isnap += 1;
    *tsnap = (*isnap)*(tend/double(snaps-1));
}

void OdeBase::before_solve (std::string dirout) {
    /* virtual, left to derived class */
    (void)dirout;
}

void OdeBase::step_extra (double t_) {
    /* virtual, left to derived class, can only be used with solve_fixed() */
    (void)t_;
}

void OdeBase::snap_extra (std::string dirout, long isnap, double tsnap) {
    /* virtual, left to derived class */
    (void)dirout; (void)isnap; (void)tsnap;
}

void OdeBase::after_solve (std::string dirout) {
    /* virtual, left to derived class */
    (void)dirout;
}

void OdeBase::check_sol_integrity () {
    for (unsigned long long i=0; i<neq; i++) {
        if ( std::isnan(sol[i]) ) {
            printf("FAILURE: Element index %llu of the solution vector is NaN\n", i);
            exit(EXIT_FAILURE);
        }
        if ( !std::isfinite(sol[i]) ) {
            printf("FAILURE: Element index %llu of the solution vector is not finite\n", i);
            exit(EXIT_FAILURE);
        }
    }
}
