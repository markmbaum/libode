//! \file ode_base.cc

#include "ode_base.h"

namespace ode {

OdeBase::OdeBase (unsigned long neq, bool need_jac) {

    //default system name
    name_ = "ode";
    //whether stuff should be printed during a solve
    quiet_ = false;
    //!whether to skip writing the solution vector to file when snapping
    silent_snap_ = false;
    //number of equations/variables in ode system
    neq_ = neq;
    //time starts at zero unless modified
    t_ = 0.0;
    //time step
    dt_ = NAN;
    //number of time steps
    nstep_ = 0;
    //number of function evaluations
    neval_ = 0;
    //interval of solution integrity checks
    icheck_ = 100;
    //array of solution values for a given time
    sol_ = new double[neq];
    //store flag
    need_jac_ = need_jac;
    //counter for jacobian evaluations
    nJac_ = 0;
    //only allocate space for Jacobian if it's needed
    if (need_jac_) {
        //allocate space for the jacobian of the system
        Jac_ = new double*[neq];
        for (unsigned long i=0; i<neq; i++) Jac_[i] = new double[neq];
        //arrays for numerical Jacobian
        f_ = new double[neq];
        g_ = new double[neq];
        //adjustment factors for numerical jacobian
        absjacdel_ = 1e-8;
        reljacdel_ = 1e-8;
    }
}

OdeBase::~OdeBase () {
    delete [] sol_;
    if (need_jac_) {
        for (unsigned long i=0; i<neq_; i++) delete [] Jac_[i];
        delete [] Jac_;
        delete [] f_;
        delete [] g_;
    }
}

//----------------------------------------------------------------------------
//use finite dif to approximate Jacobian unless virtual function reimplemented

void OdeBase::ode_jac (double *solin, double **Jout) {

    //check space was allocated for Jacobian
    if ( !need_jac_ ) ode_print_exit("the Jacobian matrix wasn't allocated for whichever solver was chosen (need_jac_ == false)");

    unsigned long i,j;
    double delsol;

    //get derivatives at current position
    ode_fun_(solin, f_);
    //perturb and calculate finite differences
    for (i=0; i<neq_; i++) {
        //perturb
        delsol = ode_max2(absjacdel_, solin[i]*reljacdel_);
        solin[i] += delsol;
        //evaluate function
        ode_fun_(solin, g_);
        //compute approximate derivatives
        for (j=0; j<neq_; j++) Jout[j][i] = (g_[j] - f_[j])/delsol;
        //unperturb
        solin[i] -= delsol;
    }
}

//-----------------------------------
//essential virtual function wrappers

void OdeBase::ode_fun_ (double *solin, double *fout) {

    //call the system of equations
    ode_fun(solin, fout);
    //increment the counter
    neval_++;
}

void OdeBase::ode_jac_ (double *solin, double **Jout) {

    //call the Jacobian function
    ode_jac(solin, Jout);
    //increment the counter
    nJac_++;
}

void OdeBase::step (double dt, bool extra) {

    //store the time step for access in derived classes
    dt_ = dt;
    //call the bare stepper
    step_(dt);
    //increment the time, a convenience variable
    t_ += dt;
    //increment the counter
    nstep_++;
    //check for nans and infs
    if (nstep_ % icheck_ == 0) check_sol_integrity();
    //do any extra stuff
    if (extra) after_step(t_);
}

//--------------
//solver support

void OdeBase::snap (std::string dirout, long isnap, double tsnap) {

    if (silent_snap_) {
        //progress report
        if (!quiet_) printf("    snap %li reached\n", isnap);
        //do any extra snapping stuff
        after_snap(dirout, isnap, tsnap);
    } else {
        //output filename
        std::string fnout = dirout + "/" + name_ + "_snap_" + ode_int_to_string(isnap);
        //write output
        ode_write(fnout.data(), sol_, neq_);
        //progress report
        if (!quiet_) printf("    snap %li written\n", isnap);
        //do any extra snapping stuff
        after_snap(dirout, isnap, tsnap);
    }
}

bool OdeBase::solve_done (double dt, double tend) {

    //if the next time step will put the time very close to or above the
    //stopping time, the main iteration is finished and a final time step
    //can be taken to get exactly to the stopping time
    if ( (t_ + dt*1.01) >= tend )
        return(true);
    return(false);
}

void OdeBase::check_sol_integrity () {

    for ( unsigned long i=0; i<neq_; i++ ) {
        //any nans?
        if ( std::isnan(sol_[i]) ) ode_print_exit("solution vector contains NaN");
        //any infs?
        if ( !std::isfinite(sol_[i]) ) ode_print_exit("solution vector contains Inf");
    }
}

void OdeBase::check_pre_solve (double tint, double dt) {

    //shouldn't be solving backward
    if (tint <= 0.0) ode_print_exit("tint must be greater than zero");
    //shouldn't be solving backward
    if (dt <= 0.0) ode_print_exit("dt must be greater than or equal to 0");
}

void OdeBase::check_pre_snaps (double dt, double *tsnap, unsigned long nsnap) {

    if (dt <= 0) ode_print_exit("dt must be greater than or equal to 0");
    if (nsnap <= 1) ode_print_exit("nsnap must be greater than 1");
    for (unsigned long i=0; i<nsnap; i++) {
        if (tsnap[i] < t_) ode_print_exit("snap times must be greater than or equal to the current time");
        if (i > 0) if (tsnap[i] <= tsnap[i-1]) ode_print_exit("snap times must monotonically increase");
    }
}

//----------------
//solver functions

void OdeBase::solve_fixed_ (double tint, double dt, bool extra) {

    //stopping time
    double tend = t_ + tint;
    //solve to within a time step of tend
    while ( !solve_done(dt, tend) ) step(dt, extra);
    //solve a fractional time step to finish
    step(tend - t_);
}

void OdeBase::solve_fixed (double tint, double dt, bool extras) {

    //checks
    check_pre_solve(tint, dt);
    if (extras) {
        //extra initial things (to be overridden in derived class)
        before_solve();
        //solve
        solve_fixed_(tint, dt, true);
        //extra completion things (to be overridden in derived class)
        after_solve();
    } else {
        solve_fixed_(tint, dt, false);
    }
}

void OdeBase::solve_fixed (double tint, double dt, const char* dirout, int inter) {

    //checks
    check_pre_solve(tint, dt);
    if (inter < 1) ode_print_exit("inter must be greater than or equal to 1");
    //indices
    unsigned long i, j;
    //output file strings
    std::string fnout;
    dirout_ = dirout;
    //stopping time
    double tend = t_ + tint;

    //output vectors
    std::vector<double> tout;
    std::vector< std::vector<double> > solout;
    solout.resize(neq_);

    //extra initial things
    before_solve();

    //store values at time zero
    for (i=0; i<neq_; i++) solout[i].push_back( sol_[i] );
    tout.push_back( t_ );
    //solve while storing solution values
    j = 0;
    while ( !solve_done(dt, tend) ) {
        //take a step
        step(dt);
        //capture values
        if ( j % inter == 0 ) {
            after_capture(t_);
            for (i=0; i<neq_; i++) solout[i].push_back( sol_[i] );
            tout.push_back( t_ );
        }
        j++;
    }
    //take a final step
    step(tend - t_);
    //store final values
    for (i=0; i<neq_; i++) solout[i].push_back( sol_[i] );
    tout.push_back( t_ );

    //write output
    for (i=0; i<neq_; i++) {
        fnout = dirout_ + "/" + name_ + "_" + ode_int_to_string(i);
        ode_write(fnout.c_str(), solout[i].data(), solout[i].size());
    }
    fnout = dirout_ + "/" + name_ + "_t";
    ode_write(fnout.c_str(), tout.data(), tout.size());

    //extra completion things (to be overridden in derived class)
    after_solve();

    //clear output directory
    dirout_ = "";
}

void OdeBase::solve_fixed (double tint, double dt, unsigned long nsnap, const char *dirout) {

    //checks
    check_pre_solve(tint, dt);
    //index
    unsigned long i;
    //stopping time
    double tend = t_ + tint;
    //make an array of the snap times
    double *tsnap = new double[nsnap];
    for (i=0; i<nsnap; i++) tsnap[i] = t_ + i*(tend - t_)/double(nsnap - 1);
    //call the other snapping function with the computed snap times
    solve_fixed(dt, tsnap, nsnap, dirout);
    //free the snap times
    delete [] tsnap;
}

void OdeBase::solve_fixed (double dt, double *tsnap, unsigned long nsnap, const char *dirout) {

    //checks
    check_pre_snaps(dt, tsnap, nsnap);
    //index
    unsigned long i;
    //store the output directory
    dirout_ = dirout;

    //extra initial things (to be overridden in derived class)
    before_solve();

    //solve to each snap time and take a snap
    for (i=0; i<nsnap; i++) {
        solve_fixed_(tsnap[i] - t_, dt);
        snap(dirout, i, t_);
    }
    //write the snap times
    if ( !silent_snap_ )
        ode_write((dirout_ + "/" + name_ + "_snap_t").data(), tsnap, nsnap);

    //extra completion things
    after_solve();

    //clear output directory
    dirout_ = "";
}

//-----
//reset

void OdeBase::reset (double t, double *sol) {

    t_ = t;
    for (unsigned long i=0; i<neq_; i++) sol_[i] = sol[i];
}

//------
//extras

void OdeBase::before_solve () { /* virtual, left to derived class */ }

void OdeBase::after_step (double t) {
    /* virtual, left to derived class, can only be used with solve_fixed() */
    (void)t;
}

void OdeBase::after_capture (double t) {
    /* virtual, left to derived class, can only be used with solve_fixed() */
    (void)t;
}

void OdeBase::after_snap (std::string dirout, long isnap, double t) {
    /* virtual, left to derived class */
    (void)dirout;
    (void)isnap;
    (void)t;
}

void OdeBase::after_solve () { /* virtual, left to derived class */ }

} // namespace ode
