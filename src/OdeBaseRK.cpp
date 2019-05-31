#include "OdeBaseRK.hpp"

OdeBaseRK::OdeBaseRK (unsigned long long neq_) : OdeBase (neq_) {

    //flag for fsal methods, assume not fsal until told otherwise
    fsal = false;
    //array of solution values for a given time
    soltemp = new double[neq];
}

//destructor
OdeBaseRK::~OdeBaseRK () {
    delete [] soltemp;
}

bool OdeBaseRK::solve_done (double dt_, double tend) {
    if ((t + dt_ >= tend) || ode_is_close(tend, t + dt_, 1e-12)) return(true);
    return(false);
}

void OdeBaseRK::solve_fixed (double tend, double dt_) {

    //store the time step
    dt = dt_;

    //get fsal methods started
    if (fsal) { ode_funk(t, sol, klast); neval++; }

    //solve to within a time step of tend
    while ( !solve_done(dt, tend) ) {
        step(dt); step_extra(t);
        //check the solution for nans and infs
        if (ncheck > 0) if (nstep % ncheck == 0) check_sol_integrity();
    }
    //solve a fractional time step to finish
    dt = tend - t;
    step(dt); step_extra(t);
}

void OdeBaseRK::solve_fixed (double tend, double dt_, const char *dirout_) {

    //index
    unsigned long long i;

    //number of steps to be taken
    unsigned long long nstep = int(ceil((tend - t)/dt_)) + 1;

    //output variables
    std::string dirout = dirout_;
    std::string fnout;
    double *tout = new double[nstep];
    double **solout = new double*[neq];
    for (i=0; i<neq; i++) solout[i] = new double[nstep];

    //store the time step
    dt = dt_;

    //extra initial things (to be overridden in derived class)
    before_solve(dirout);

    //get fsal methods started
    if (fsal) { ode_funk(t, sol, klast); neval++; }

    //store values at time zero
    for (i=0; i<neq; i++) { solout[i][0] = sol[i]; } tout[0] = t;
    //solve while storing every solution value
    unsigned long long count = 1;
    while ( !solve_done(dt, tend) ) {
        step(dt); step_extra(t);
        //store values
        for (i=0; i<neq; i++) { solout[i][count] = sol[i]; } tout[count] = t;
        //check the solution for nans and infs
        if (ncheck > 0) if (nstep % ncheck == 0) check_sol_integrity();
        //update the counter
        count++;
    }
    //take a final step to the finish line
    dt = tend - t;
    step(dt); step_extra(t);

    //store final values
    for (i=0; i<neq; i++) { solout[i][count] = sol[i]; } tout[count] = t;

    //write output
    for (i=0; i<neq; i++) {
        fnout = dirout + "/sol_" + std::to_string(i);
        ode_write_double(fnout.data(), solout[i], nstep);
    }
    fnout = dirout + "/t";
    ode_write_double(fnout.data(), tout, nstep);

    //extra completion things (to be overridden in derived class)
    after_solve(dirout);

    delete [] tout;
    for (i=0; i<neq; i++) { delete [] solout[i]; } delete [] solout;
}

void OdeBaseRK::solve_fixed (double tend, double dt_, int snaps, const char *dirout_) {

    //snapping variables
    std::string dirout = dirout_;
    std::vector<double> tout;
    long isnap = 0;
    double tsnap = 0.0;

    //store the time step
    dt = dt_;

    //extra initial things (to be overridden in derived class)
    before_solve(dirout);

    //get fsal methods started
    if (fsal) { ode_funk(t, sol, klast); neval++; }

    //write snapshot at time zero
    snap(dirout, &isnap, &tsnap, tend, snaps, &tout);
    //solve while writing snaps
    while ( !( solve_done(dt, tend) && ode_is_close(tsnap, tend, 1e-12) ) ) {
        //if the snap time is before the next time step, step to the snap time
        if ( ode_is_close(tsnap, t + dt, 1e-12) || ( tsnap < t + dt ) ) {
            dt = tsnap - t;
            step(dt); step_extra(t);
            snap(dirout, &isnap, &tsnap, tend, snaps, &tout);
            dt = dt_;
        } else {
            step(dt); step_extra(t);
        }
        //check the solution for nans and infs
        if (ncheck > 0) if (nstep % ncheck == 0) check_sol_integrity();
    }
    dt = tend - t;
    step(dt); step_extra(t);
    snap(dirout, &isnap, &tsnap, tend, snaps, &tout);
    ode_write_double((dirout + "/snap_times").data(), tout.data(), snaps);

    //extra completion things (to be overridden in derived class)
    after_solve(dirout);
}
