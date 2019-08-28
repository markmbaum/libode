//! \file ode_base.h
/*!
\mainpage libode

libode is a library of C++ classes for solving systems of ordinary differential equations in autonomous form. All of the solvers are single-step, Runge-Kutta-like methods. There are explicit, adaptive solvers up to the ninth order. The repository also includes Rosenbrock methods, a singly-diagonal implicit Runge-Kutta (SDIRK) method, and several fully implicit Runge-Kutta methods. With the current collection of solvers and features, `libode` is well suited to any non-stiff systems and to stiff systems that are tightly coupled and have a known Jacobian (ones that don't require sparse or banded matrix routines).

Several of the solvers and much more detail on the methods can be found in these amazing books:
+ Hairer, E., Nørsett, S. P. & Wanner, G. Solving Ordinary Differential Equations I: Nonstiff Problems. (Springer-Verlag, 1987).
+ Hairer, E. & Wanner, G. Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems. (Springer, 1996).

\section sec_compiling Compiling

`libode` was written to provide easy access to class-based ODE solvers without dependencies or specialized compiling processes. There is only one step to take before compiling. Consequently, the library is slim on features and doesn't provide access to things like sparse matrices. For many systems of ODEs, though, `libode` should make it easy to build an integrator and enjoy the speed of C++ and [openmp](https://en.wikipedia.org/wiki/OpenMP) without the headaches of large, complex packages.

First, before any of the `libode` classes can be compiled, you must copy the `_config.mk` file to `config.mk` and edit that file to specify the compiler settings you'd like the Makefile to use. This shouldn't be complicated. If you are using a current version of the GNU C++ compiler (g++), the contents of the template config file can likely be used without modification. There are also commented lines for use with the Intel C++ compiler (icpc), if that is available. To compile all the classes, simply run `make` in the top directory.

The Makefile compiles all of the necessary code into the `obj` folder, then archives it in the `bin` directory as a file called `libode.a`. To use the solvers, you can link `libode.a` (in the `bin` directory) or the object files directly (in the `obj` directory) when compiling your derived class. You must also the the header files in the `src` directory, as there is not a single header file for the library. All of the classes have their header file name displayed in the documentation. Linking the solver classes requires something like

`-I<path>/libode/src -L<path>/libode/bin -lode`

when compiling derived code, with `<path>` replaced by path elements leading to the libode directory. For some examples of how to link a derived class to `libode` and create a program to run integrations, see the examples folder.

Test programs are compiled with `make tests` and they can all be run in sequence with the `run_all_tests.sh` script (which uses Python to plot the test results).

\section sec_usage Using the Solvers

To integrate a specific system of ODEs, a new class must be created to inherit from one of the solver classes. This new inheriting class must
1. Define the system of ODEs to be solved by implementing the `ode_fun()` function. This is a virtual function in the base classes. Once it is implemented, it can be used by the stepping and solving functions.
2. Set initial conditions using the `set_sol()` function.
3. Optionally implement the `ode_jac()` function for implicit methods. This is also a virtual function in the base classes. If it's not overridden, a finite-difference estimate of the Jacobian is used.

For flexibility, the derived class could be a template, so that the solver/method can be chosen when the class is constructed. Other than defining the system of equations and setting initial conditions, the derived class can store whatever information and implement whatever other methods are necessary. This could be something simple like an extra function for setting initial conditions. It could, however, comprise any other system that needs to run on top of an ODE solver, like the spatial discretization of a big PDE solver.

Each solver has a `step` method that can be used to integrate a single step with a specified step size. Each solver class also has a `solve_fixed()` method and, if it's an adaptive class, a `solve_adaptive()` method. These functions return nothing and both have the same four call signatures:

1. `void solve_fixed (double tint, double dt)`

   Simply advances the solution for a specified length of the independent variable. The independent variable is assumed to be time, so `tint` is the integration time and `dt` is the time step to use (or the initial time step for adaptive solves).

2. `void solve_fixed (double tint, double dt, const char *dirout, int inter)`

   Integrates for a duration of `tint` using time step (or initial time step) `dt` and writes solution values after every `inter` steps to the directory `dirout`. For example, if `inter` is one, the solution at every step is written to file. If `inter` is two, every other step is written.

3. `void solve_fixed (double tint, double dt, unsigned long nsnap, const char *dirout)`

   Integrates and writes `nsnap` even spaced snapshots of the solution into the directory `dirout`.

4. `void solve_fixed (double dt, double *tsnap, unsigned long nsnap, const char *dirout)`

   Integrates and writes snapshots at the times specified in `tsnap` into the directory `dirout`.
*/

#ifndef ODE_BASE_H_
#define ODE_BASE_H_

#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "ode_io.h"
#include "ode_util.h"

//!Lowest base class for all solvers
/*!
This is the deepest base class, upon which all integrators and other base classes are built. It provides basic variables like the solution array, the independent variable, and a string defining the integrator's name. It implements the solve_fixed function for integrating with a fixed time step. When constructing, it will allocate space for a Jacobian matrix if the need_jac flag is true.
*/
class OdeBase {

    public:
        //!constructs
        /*!
        \param[in] neq number of equations in ODE system
        \param[in] need_jac flag signaling whether the Jacobian of the system is needed
        */
        OdeBase (unsigned long neq, bool need_jac);
        //!destructs
        ~OdeBase ();

        //-------------------
        //getters and setters

        //!gets the name of the ODE system
        const char *get_name () { return(name_.c_str()); }
        //!gets the name of the solver/method
        const char *get_method () { return(method_.c_str()); }
        //!gets the size of the ODE system
        unsigned long get_neq () { return(neq_); }
        //!gets the current value of the independent variable
        double get_t () { return(t_); }
        //!gets the most recent time step size
        double get_dt () { return(dt_); }
        //!gets the solution array
        double *get_sol () { return(sol_); }
        //!gets an element of the solution array
        double get_sol (unsigned long i) { return(sol_[i]); }
        //!gets the total number of steps taken
        unsigned long long get_nstep () { return(nstep_); }
        //!gets the total number of ODE system evaluation
        unsigned long long get_neval () { return(neval_); }
        //!gets the number of steps after which the solution is checked for integrity
        unsigned long long get_icheck () { return(icheck_); }
        //!gets the total number of Jacobian evaluations performed
        unsigned long long get_nJac () { return(nJac_); }

        //!sets an element of the solution array
        void set_sol (unsigned long i, double x) { sol_[i] = x; }
        //!copies an array into the solution array
        void set_sol (double *sol) { for(unsigned long i=0; i<neq_; i++) sol_[i] = sol[i]; }
        //!sets the name of the ODE system
        void set_name (std::string name) { name_ = name; }
        //!sets the name of the ODE system
        void set_name (const char *name) { name_ = name; }
        //!sets the number of steps after which the solution is checked for integrity
        void set_icheck (unsigned long icheck) { icheck_ = icheck; }

        //----------------
        //solver functions

        //!increments the step counter and the time, checks the solution integrity if needed, stores the time step in the object, and executes after_step()
        void step (double dt);

        //!integrates for a specified duration of independent variable without output
        /*!
        \param[in] tint total integration time
        \param[in] dt time step size
        */
        void solve_fixed (double tint, double dt);

        //!lots of output, solves and stores every "inter" point along the way
        /*!
        \param[in] tint total integration time
        \param[in] dt time step size
        \param[in] dirout output directory (must already exist)
        \param[in] inter interval of steps to store and output
        */
        void solve_fixed (double tint, double dt, const char *dirout, int inter);

        //!solves and writes evenly spaced snapshots
        /*!
        \param[in] tint total integration time
        \param[in] dt time step size
        \param[in] nsnap number of snapshots to output
        \param[in] dirout output directory (must already exist)
        */
        void solve_fixed (double tint, double dt, unsigned long nsnap, const char *dirout);

        //!solves and writes snapshots at times specified in the tsnap array
        /*!
        \param[in] dt time step size
        \param[in] tsnap array of desired snapshot times
        \param[in] nsnap number of snapshots (length of tsnap)
        \param[in] dirout output directory (must already exist)
        */
        void solve_fixed (double dt, double *tsnap, unsigned long nsnap, const char *dirout);

        //!reset to a specified time and initial condition array
        /*!
        \param[in] t new value of independent variable
        \param[in] sol array of new values for solution array
        */
        void reset (double t, double *sol);

    protected:

        //!integrates without output or any counters, trackers, extra functions...
        /*!
        \param[in] tint total integration time
        \param[in] dt time step size
        */
        void solve_fixed_ (double tint, double dt);

        //----------------------
        //basic solver variables

        //!the "name" of the system, which is used for output
        std::string name_;
        //!the "name" of the solver/method, as in "Euler" or "RK4"
        std::string method_;
        //!number of equations in the system of ODEs
        unsigned long neq_;
        //!time, initialized to zero
        double t_;
        //!time step is stored and updated during solves
        double dt_;
        //!array for the solution, changing over time
        double *sol_;
        //!number of time steps
        unsigned long long nstep_;
        //!function evaluation counter, must be incremented in step() when defined
        unsigned long long neval_;
        //!interval of steps after which to check for nans and infs (zero to ignore)
        unsigned long long icheck_;
        //!storage for the ODE system's Jacobian matrix
        double **Jac_;
        //!counter for jacobian evaluations
        unsigned long long nJac_;
        //!absolute adjustment fraction for numerical Jacobian, if needed
        double absjacdel_;
        //!relative adjustment fraction for numerical Jacobian,`` if needed
        double reljacdel_;

        //---------------------------
        //essential virtual functions

        //!evaluates the system of ODEs in autonomous form and must be defined by a derived class
        /*!
        \param[in] solin current solution array
        \param[out] fout evaluation of system of ordinary differential equations
        */
        virtual void ode_fun (double *solin, double *fout) = 0;

        //!evaluates the system's Jacobian matrix, also in autonomous form, and can either be defined in a derived class or left to numerical approximation
        /*!
        \param[in] solin current solution array
        \param[out] Jout Jacobian of ode_fun
        */
        virtual void ode_jac (double *solin, double **Jout);

        //!advances a single time step (without changing counters or the time) and must be defined in the derived class implementing the solver/method
        /*!
        \param[in] dt time step size
        */
        virtual void step_ (double dt) = 0;

        //---------------------------------------
        //wrappers of essential virtual functions

        //!wrapper, calls ode_fun() and increases the neval counter by one
        void ode_fun_ (double *solin, double *fout);
        //!wrapper, calls ode_jac() and increments nJac;
        void ode_jac_ (double *solin, double **Jout);

        //------
        //extras

        //!does any extra stuff before starting a solve
        virtual void before_solve ();
        //!does any extra stuff after each step
        /*!
        \param[in] t current value of ODE system's independent variable
        */
        virtual void after_step (double t);
        //!does any extra stuff only when a step is captured
        /*!
        \param[in] t current value of ODE system's independent variable
        */
        virtual void after_capture (double t);
        //!does any extra stuff after each snap
        /*!
        \param[in] dirout output directory (must already exist)
        \param[in] isnap index of snap being taken (from 0)
        \param[in] t current value of ODE system's independent variable
        */
        virtual void after_snap (std::string dirout, long isnap, double t);

        //!does any extra stuff after completing a solve
        virtual void after_solve ();

        //--------------
        //solver support

        //!writes the current value of the solution to a binary file
        /*!
        \param[in] dirout output directory (must already exist)
        \param[in] isnap the index of the snap
        \param[in] tsnap independent variable value for snapshot
        */
        void snap (std::string dirout, long isnap, double tsnap);

        //!checks if the solution is within a single time step of the end point
        /*!
        \param[in] dt time step duration
        \param[in] tend end time
        \return true if the solution is within a step of the end, otherwise false
        */
        bool solve_done (double dt, double tend);

        //!checks solution for nans and infs, exiting the program if they're found
        void check_sol_integrity ();

        //!checks that a solve can be performed with given tend and dt values
        /*
        \param[in] tint duration of integration
        \param[in] dt time step size
        */
        void check_pre_solve (double tint, double dt);

        //!checks that snap times are monotonically increasing and > current time
        /*
        \param[in] dt time step size
        \param[in] tsnap array of snapshot times during integration
        \param[in] nsnap number of snapshots (length of tsnap)
        */
        void check_pre_snaps (double dt, double *tsnap, unsigned long nsnap);

    private:

        //flag for whether the Jacobian is being used
        bool need_jac_;
        //arrays for evaluating numerical jacobian
        double *f_, *g_, *soljac_;
};

#endif
