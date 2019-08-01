/*
This is the deepest base class, upon which all integrators and other base classes are built. It provides basic variables like the solution array, the independent variable, and a string defining the integrator's name. It implements the solve_fixed function for integrating with a fixed time step. When constructing, it will allocate space for a Jacobian matrix if the need_jac flag is true.
*/

#ifndef ODE_BASE_H_
#define ODE_BASE_H_

#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "ode_io.h"
#include "ode_util.h"

class OdeBase {

    public:
        //constructor
        OdeBase (unsigned long neq, bool need_jac);
        //destructor
        ~OdeBase ();

        //-------------------
        //getters and setters

        const char *get_name () { return(name_.c_str()); }
        const char *get_method () { return(method_.c_str()); }
        unsigned long get_neq () { return(neq_); }
        double get_t () { return(t_); }
        double get_dt () { return(dt_); }
        double *get_sol () { return(sol_); }
        double get_sol (unsigned long i) { return(sol_[i]); }
        unsigned long long get_nstep () { return(nstep_); }
        unsigned long long get_neval () { return(neval_); }
        unsigned long long get_ncheck () { return(ncheck_); }
        unsigned long long get_nJac () { return(nJac_); }

        void set_sol (unsigned long i, double z) { sol_[i] = z; }
        void set_sol (double *sol) { for(unsigned long i=0; i<neq_; i++) sol_[i] = sol[i]; }
        void set_name (std::string name) { name_ = name; }
        void set_name (const char *name) { name_ = name; }
        void set_ncheck (unsigned long ncheck) { ncheck_ = ncheck; }

        //increments the step counter and the time, checks the solution integrity,
        //stores the time step in the object, and executes after_step()
        void step (double dt);

        //----------------
        //solver functions

        //solves to a particular time with a fixed time step
        //    tend - stopping time, time at end of solution
        //    dt - time step
        //    dirout_ - relative path to preexisting directory for output
        //    inter - interval of steps to store, i.e. if s=3 every third step is stored
        //    nsnap - number of snapshots to output
        //    tsnap - array of specific times to take snapshots

        //no output, solves to tend and executes extras like before_solve(), ...
        void solve_fixed (double tint, double dt);
        //lots of output, solves and stores every "inter" point along the way
        void solve_fixed (double tint, double dt, const char *dirout, int inter);
        //some output, solves and writes evenly spaced snapshots
        void solve_fixed (double tint, double dt, unsigned long nsnap, const char *dirout);
        //some output, solves and writes snapshots at times in the tsnap array
        void solve_fixed (double dt, double *tsnap, unsigned long nsnap, const char *dirout);

        //reset to a certain time and solution array
        void reset (double t, double *sol);

    protected:
        //no output, solves to tend without before_solve(), after_solve(), ...
        void solve_fixed_ (double tint, double dt);

        //----------------------
        //basic solver variables

        //the "name" of the system, which is used for output
        std::string name_;
        //the "name" of the solver/method, as in "Euler" or "RK4"
        std::string method_;
        //number of equations in the system of ODEs
        unsigned long neq_;
        //time, initialized to zero
        double t_;
        //time step is stored and updated during solves
        double dt_;
        //array for the solution, changing over time
        double *sol_;
        //number of time steps
        unsigned long long nstep_;
        //function evaluation counter, must be incremented in step() when defined
        unsigned long long neval_;
        //interval of steps after which to check for nans and infs (zero to ignore)
        unsigned long long ncheck_;
        //storage for the ODE system's Jacobian matrix
        double **Jac_;
        //counter for jacobian evaluations
        unsigned long long nJac_;
        //adjustment fractions for numerical jacobian
        double absjacdel_, reljacdel_;

        //---------------------------
        //essential virtual functions

        /*evaluates the system of ODEs, in AUTONOMOUS form, and must be \
        defined by a derived class
        input:
            solin - current solution array
        output:
            fout - evaluation of system of ordinary differential equations*/
        virtual void ode_fun (double *solin, double *fout) = 0;

        /*evaluates the system's Jacobian matrix, also in autonomous form, and
        can either be defined in a derived class or left to numerical
        approximation
        input:
            solin - current solution array
        output:
            Jout - finite dif approximation to Jacobian of ode_fun*/
        virtual void ode_jac (double *solin, double **Jout);

        /*advances a single time step (without changing counters or the time) and
        must be defined in the derived class implementing the solver/method
        input:
            dt - time step size*/
        virtual void step_ (double dt) = 0;

        //---------------------------------------
        //wrappers of essential virtual functions

        //also increases the neval counter by one
        void ode_fun_ (double *solin, double *fout);

        //wrapper to increment nJac;
        void ode_jac_ (double *solin, double **Jout);

        //------
        //extras

        //does any extra stuff before starting a solve
        virtual void before_solve ();
        //does any extra stuff after each step
        virtual void after_step (double t);
        //does any extra stuff only when a step is captured
        virtual void after_capture (double t);
        //does any extra stuff after each snap
        virtual void after_snap (std::string dirout, long isnap, double t);
        //does any extra stuff after completing a solve
        virtual void after_solve ();

        //--------------
        //solver support

        /*writes the current value of the solution to a binary file
        input:
            dirout - directory to write solution array into
            isnap - the index of the snap
            tsnap - the current time or time of snapshot*/
        void snap (std::string dirout, long isnap, double tsnap);

        /*checks if the solution is within a single time step of the end point
        input:
            dt_ - time step duration
            tend - end time
        output:
            true if the solution is within a step of the end, otherwise false*/
        bool solve_done (double dt, double tend);

        //checks solution for nans and infs, exiting the program if they're found
        void check_sol_integrity ();

        //checks that a solve can be performed with given tend and dt values
        void check_pre_solve (double tint, double dt);

        //checks that snap times are monotonically increasing and > current time
        void check_pre_snaps (double dt, double *tsnap, unsigned long nsnap);

    private:
        //flag for whether the Jacobian is being used
        bool need_jac_;
        //arrays for evaluating numerical jacobian
        double *f_, *g_, *soljac_;
};

#endif
