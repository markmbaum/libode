/*
The OdeAdaptive class is a base for any integrator that can modify its own time step as it integrates and reject steps. The class inherits from the deepest base class, ode_base. It implements three virtual functions to enable time step adaptation:
  1) The 'adapt' function is meant to execute whatever calculations are needed
     to subsequently assess whether the current time step should be rejected
     and what the next time step should be.
  2) The 'is_rejected' function simply furnishes a boolean answering the
     question, "should the step just taken be rejected?"
  3) The 'dt_next' function furnishes the next time step, whether the step
     current step was rejected or not.
These functions are used with local error tolerances to integrate with an
efficiently adaptive time step. Specifically, they are used in the step_adaptive_ method of this class, which is called in each of the solve_adaptive methods.
*/

#ifndef ODE_ADAPTIVE_H_
#define ODE_ADAPTIVE_H_

#include <cstdio>
#include <cmath>
#include <vector>
#include <cstring>
#include <string>

#include "ode_io.h"
#include "ode_util.h"

#include "ode_base.h"

class OdeAdaptive : public OdeBase {

    public:
        //constructor
        OdeAdaptive (unsigned long neq, bool need_jac);
        //destructor
        ~OdeAdaptive ();

        //-------------------
        //getters and setters

        unsigned long long get_nrej () { return(nrej_); }
        double get_abstol () { return(abstol_); }
        double get_reltol () { return(reltol_); }
        double get_dtmax () { return(dtmax_); }

        void set_abstol (double tol) { abstol_ = tol; };
        void set_reltol (double tol) { reltol_ = tol; };
        void set_tol (double tol) { abstol_ = tol; reltol_ = tol; };
        void set_dtmax (double dtmax) { dtmax_ = dtmax; }

        //integrates with adaptive time stepping
        //    tend - stopping time, time at end of solution
        //    dt - time step
        //    snaps - number of snapshots to output
        //    dirout - relative path to preexisting directory for output

        //no output, solves to tend and executes extras like before_solve(), ...
        void solve_adaptive (double tint, double dt0);
        //lots of output, solves and stores every "inter" point along the way
        void solve_adaptive (double tint, double dt0, const char *dirout, int inter);
        //some output, solves and writes evenly spaced snapshots
        void solve_adaptive (double tint, double dt0, unsigned long nsnap, const char *dirout);
        //some output, solves and writes snapshots at times in the tsnap array
        void solve_adaptive (double dt0, double *tsnap, unsigned long nsnap, const char *dirout);

    protected:
        //no output, solves to tend without before_solve(), after_solve(), ...
        void solve_adaptive_ (double tint, double dt0);

        //------------------------
        //basic adaptive variables

        //counter for rejected steps
        unsigned long long nrej_;
        //absolute and relative error tolerances
        double abstol_, reltol_;
        //maximum allowable time step
        double dtmax_;

        //executes whatever calculations need to be performed for adapting
        //this should include a determination of whether a step is rejected
        //and a calculation of the next time step size
        virtual void adapt (double abstol, double reltol) = 0;
        //retreives a bool determining whether a step is accepted/rejected
        virtual bool is_rejected () = 0;
        //retrieves the best time step for the next step
        virtual double dt_next () = 0;

        //determines whether an adaptive solve is finished
        bool solve_done_adaptive (double tend);

        //wrappers
        bool step_adaptive_ (double dt);
        double dt_next_ (double tend);

    private:
        //previous solution, in case a step is rejected
        double *solprev_;
};

#endif
