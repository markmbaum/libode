#include "OdeBaseARK.hpp"

OdeBaseARK::OdeBaseARK (long int neq_, int errord_) : OdeBaseRK (neq_) {

    //order of the LOWER order solution used for error estimation
    errord = errord_;

    //low order solution for error estimation
    sollo = new double[neq];
    //high order solution for error estimation
    solhi = new double[neq];
    //storage of previous y in case of rejection
    solprev = new double[neq];
    //variables for adaptive time stepping and step rejection
    dtmax = std::numeric_limits<double>::infinity(); //maximum time step size
    facsafe = 0.9;      //fraction to multiply facopt by to be cautious
    facmin = 1.0/100.0; //minumum time step reduction factor
    facmax = 10.0;      //maximum time step increase factor
    abstol = 1e-6;      //default absolute error tolerance
    reltol = 1e-6;      //default relative error tolerance
    nrej = 0;           //number of rejected steps
}

//destructor
OdeBaseARK::~OdeBaseARK () {
    delete [] sollo;
    delete [] solhi;
    delete [] solprev;
}

double OdeBaseARK::error (double abstol, double reltol) {

    //index
    unsigned long long i;
    //storage for evaluating error terms
    double sc, d, err=0.0;

    //compute the 2 norm of the difference between higher and lower order steps
    //with a scaling factor that combines relative and absolute error
    for (i=0; i<neq; i++) {
        d = fabs(sollo[i] - solhi[i]);
        sc = abstol + reltol*ode_max2(solprev[i], sol[i]);
        err += (d/sc)*(d/sc);
    }
    err = sqrt(err/neq);
    return(err);
}

double OdeBaseARK::facopt (double err) {
    //compute the "optimal" adjustment factor for the time step
    return(
        ode_min2(facmax, ode_max2(facmin, facsafe*pow(1.0/err, 1.0/(errord + 1.0))))
    );
}

void OdeBaseARK::reverse_step () {
    //put the previous step's solution back into the current solution array
    for (unsigned long long i=0; i<neq; i++) sol[i] = solprev[i];
    //reverse the time and the counters
    t -= dt; nrej++; nstep--;
    //if fsal, recompute the final k values, which become the k1 values
    if (fsal) { ode_funk(t, sol, klast); neval++; }
}

void OdeBaseARK::solve_adaptive (double tend, double dt0) {

    //error storage for adapting
    double err;
    //store the initial time step
    dt = dt0;

    //if the method has the FSAL property, initialize k5 as if it were k1
    if (fsal) { ode_funk(t, sol, klast); neval++; }

    //solve to tend
    while ( !ode_is_close(t, tend, 1e-8) && (t < tend) ) {
        //take a step
        step(dt);
        //compute error estimate
        err = error(abstol, reltol);
        //if error is too large, reverse the step to try again
        if (err > 1.0) reverse_step();
        //new step size shouldn't put t over tend
        dt = ode_min2(dt*facopt(err), tend - t);
        //new step size also shouldn't exceed dtmax
        dt = ode_min2(dt, dtmax);
        //check the solution for nans and infs
        if (ncheck > 0) if (nstep % ncheck == 0) check_sol_integrity();
    }
}

void OdeBaseARK::solve_adaptive (double tend, double dt0, const char *dirout_) {

    //index
    unsigned long long i;
    //output variables
    std::vector<double> *solout = new std::vector<double>[neq];
    std::vector<double> tout;
    std::string dirout = dirout_;
    std::string fnout;
    //other variables
    double err;

    //store the initial time step
    dt = dt0;

    //extra initial things (to be overridden in derived class)
    before_solve(dirout);

    //if the method has the FSAL property, initialize k5 as if it were k1
    if (fsal) { ode_funk(t, sol, klast); neval++; }

    //store values at time zero
    for (i=0; i<neq; i++) solout[i].push_back(sol[i]);
    tout.push_back(t);
    //solve to the end of time (tend)
    while ( !ode_is_close(t, tend, 1e-8) && (t < tend) ) {
        //take a step
        step(dt);
        //compute error estimate
        err = error(abstol, reltol);
        //if error is too large, reverse the step to try again, otherwise store values
        if (err > 1.0) {
            reverse_step();
        } else {
            //store values
            for (i=0; i<neq; i++) solout[i].push_back(sol[i]);
            tout.push_back(t);
        }
        //new step size shouldn't put t over tend
        dt = ode_min2(dt*facopt(err), tend - t);
        //new step size also shouldn't exceed dtmax
        dt = ode_min2(dt, dtmax);
        //check the solution for nans and infs
        if (ncheck > 0) if (nstep % ncheck == 0) check_sol_integrity();
    }

    //write output
    for (i=0; i<neq; i++) {
        fnout = dirout + "/sol" + std::to_string(i);
        ode_write_double(fnout.data(), solout[i].data(), int(solout[i].size()));
    }
    fnout = dirout + "/t";
    ode_write_double(fnout.data(), tout.data(), int(tout.size()));

    //extra completion things (to be overridden in derived class)
    after_solve(dirout);
}

void OdeBaseARK::solve_adaptive (double tend, double dt0, int snaps, const char *dirout_) {

    //snapping variables
    std::string dirout = dirout_;
    std::vector<double> tout;
    double tsnap = 0.0;
    long isnap = 0;
    //other variables
    double err;

    //store the initial time step
    dt = dt0;

    //extra initial things (to be overridden in derived class)
    before_solve(dirout);

    //if the method has the FSAL property, initialize k5 as if it were k1
    if (fsal) { ode_funk(t, sol, klast); neval++; }

    //write snapshot at time zero
    snap(dirout, &isnap, &tsnap, tend, snaps, &tout);
    //solve to the end of time (tend)
    dt = (dt0 < tsnap - t) ? dt0 : tsnap - t;
    while ( !ode_is_close(t, tend, 1e-8) && (t < tend) ) {
        //take a step
        step(dt);
        //compute error estimate
        err = error(abstol, reltol);
        //if error is too large, reverse the step to try again
        if (err > 1.0) {
            reverse_step();
        } else if ( t == tsnap ) {
            snap(dirout, &isnap, &tsnap, tend, snaps, &tout);
        }
        //new step size shouldn't put t over tend
        dt = ode_min2(dt*facopt(err), tend - t);
        //new step size also shouldn't exceed dtmax
        dt = ode_min2(dt, dtmax);
        //new step shouldn't go past next snap time
        dt = ode_min2(dt, tsnap - t);
        //check the solution for nans and infs
        if (ncheck > 0) if (nstep % ncheck == 0) check_sol_integrity();
    }
    //write the times
    ode_write_double((dirout + "/snap_times").data(), tout.data(), int(tout.size()));

    //extra completion things (to be overridden in derived class)
    after_solve(dirout);
}
