//! \file ode_adaptive.cc

#include "ode_adaptive.h"

OdeAdaptive::OdeAdaptive (unsigned long neq, bool need_jac) :
    OdeBase (neq, need_jac) {

    //rejected step counter
    nrej_ = 0;
    //error tolerances
    abstol_ = 1e-6;
    reltol_ = 1e-6;
    //hard upper time step limit
    dtmax_ = INFINITY;
    //previous solution for step reversals
    solprev_ = new double[neq];
}

OdeAdaptive::~OdeAdaptive () {
    delete [] solprev_;
}

//---------------------------------------
//integration with adaptive time stepping

void OdeAdaptive::solve_adaptive_ (double tint, double dt0, bool extra) {

    //store the time step
    double dt = dt0;
    //stopping time
    double tend = t_ + tint;
    //don't step over the edge
    if (t_ + dt > tend) dt = tend - t_;
    //solve to the stopping time
    while ( !solve_done_adaptive(tend) ) {
        //take a step, which might be rejected
        step_adaptive_(dt, extra);
        //compute a new time step
        dt = dt_adapt_(tend);
    }
}

void OdeAdaptive::solve_adaptive (double tint, double dt0, bool extras) {

    //checks
    check_pre_solve(tint, dt0);
    if (extras) {
        //extra initial things
        before_solve();
        //solve
        solve_adaptive_(tint, dt0, true);
        //extra completion things (to be overridden in derived class)
        after_solve();
    } else {
        solve_adaptive_(tint, dt0, false);
    }
}

void OdeAdaptive::solve_adaptive (double tint, double dt0, const char *dirout, int inter) {

    //checks
    check_pre_solve(tint, dt0);
    if (inter < 1) ode_print_exit("inter must be greater than or equal to 1");
    //indices
    unsigned long i, j;
    //output file strings
    std::string fnout;
    dirout_ = dirout;
    //output vectors
    std::vector<double> tout;
    std::vector<double> *solout = new std::vector<double>[neq_];
    //stopping time
    double tend = t_ + tint;

    //extra initial things (to be overridden in derived class)
    before_solve();

    //store values at time zero
    for (i=0; i<neq_; i++) solout[i].push_back( sol_[i] );
    tout.push_back( t_ );
    //solve while storing solution values
    double dt = dt0;
    bool suc;
    j = 0;
    while ( !solve_done_adaptive(tend) ) {
        //take a step
        suc = step_adaptive_(dt);
        //compute a new time step
        dt = dt_adapt_(tend);
        //if the step was successful, see if values should be stored
        if (suc) {
            //count another successful step
            j++;
            //check if storing
            if ( j % inter == 0 ) {
                after_capture(t_);
                //store
                for (i=0; i<neq_; i++) solout[i].push_back( sol_[i] );
                tout.push_back( t_ );
            }
        }
    }
    //store the final time step if need be
    if ( j-1 % inter != 0 ) {
        for (i=0; i<neq_; i++) solout[i].push_back( sol_[i] );
        tout.push_back( t_ );
    }

    //write output
    for (i=0; i<neq_; i++) {
        fnout = dirout_ + "/" + name_ + "_" + ode_int_to_string(i);
        ode_write(fnout.data(), solout[i].data(), solout[i].size());
    }
    fnout = dirout_ + "/" + name_ + "_t";
    ode_write(fnout.data(), tout.data(), tout.size());

    //extra completion things (to be overridden in derived class)
    after_solve();

    //clear output directory
    dirout_ = "";
}

void OdeAdaptive::solve_adaptive (double tint, double dt0, unsigned long nsnap, const char *dirout) {
    //checks
    check_pre_solve(tint, dt0);
    //stopping time
    double tend = t_ + tint;
    //make an array of the snap times
    double *tsnap = new double[nsnap];
    for (unsigned long i=0; i<nsnap; i++) tsnap[i] = t_ + i*(tend - t_)/double(nsnap - 1);
    //call the other snapping function with the computed snap times
    solve_adaptive(dt0, tsnap, nsnap, dirout);
    //free the snap times
    delete [] tsnap;
}

void OdeAdaptive::solve_adaptive (double dt0, double *tsnap, unsigned long nsnap, const char *dirout) {

    //checks
    check_pre_snaps(dt0, tsnap, nsnap);
    //index
    unsigned long i;
    //store the output directory
    dirout_ = dirout;

    //extra initial things (to be overridden in derived class)
    before_solve();

    //solve to each snap time and take a snap
    double dt = dt0;
    for (i=0; i<nsnap; i++) {
        //basic adaptive solve
        solve_adaptive_(tsnap[i] - t_, dt);
        //write the current solution to file
        snap(dirout, i, t_);
        //store the estimate for the current best time step
        dt = dt_adapt_(INFINITY);
    }
    //write the snap times
    if ( !silent_snap_ )
        ode_write((dirout_ + "/" + name_ + "_snap_t").data(), tsnap, nsnap);

    //extra completion things (to be overridden in derived class)
    after_solve();

    //clear output directory
    dirout_ = "";
}

//--------------
//adapting stuff

void OdeAdaptive::adapt (double abstol, double reltol) {
    (void)abstol;
    (void)reltol;
}

bool OdeAdaptive::is_rejected () { return(false); }

double OdeAdaptive::dt_adapt () {
    ode_print_exit("An adaptive solving method was called without a time step choosing algorithm implemented! You must at least implement the dt_adapt() function to use solve_adaptive().");
    return(1);
}

bool OdeAdaptive::solve_done_adaptive (double tend) {

    if ( ode_is_close(t_, tend, 1e-13) || (t_ >= tend) )
        return(true);
    return(false);
}

//--------------
//wrappers

bool OdeAdaptive::step_adaptive_ (double dt, bool extra) {

    //store the current solution in case the step is rejected
    memcpy(solprev_, sol_, neq_*sizeof(double));
    //store the time step for access in derived classes
    dt_ = dt;
    //call the basic stepper
    step_(dt);
    //count the step, even if it's subsequently rejected
    nstep_++;
    //execute adaptive calculations
    adapt(abstol_, reltol_);
    //see if the step should be rejected
    if ( is_rejected() ) {
        //reject the step by putting solprev back into sol
        memcpy(sol_, solprev_, neq_*sizeof(double));
        //count the rejection
        nrej_++;
        //return failure
        return(false);
    } else {
        //increment the time
        t_ += dt;
        //check for nans and infs
        if (nstep_ % icheck_ == 0) check_sol_integrity();
        //do any extra stuff
        if (extra) after_step(t_);
        //return success
        return(true);
    }
}

double OdeAdaptive::dt_adapt_ (double tend) {

    //get next time step from virtual function
    double dt = dt_adapt();
    //make sure the new dt doesn't exceed the stopping time
    if ( tend < t_ + dt*1.01 ) dt = tend - t_;
    //make sure the new dt is not larger than the maximum time step
    if ( dt > dtmax_ ) dt = dtmax_;

    return(dt);
}
