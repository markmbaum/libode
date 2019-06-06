/*
OdeBaseARK serves as a base for single step, explicit Runge-Kutta methods with embedded solutions and adaptive time stepping (ark = adaptive explicit runge kutta). It contains variables and functions that are naturally shared by these methods, like arrays for high and low order solutions that are used to estimate error and adjust the time step. It inherits functions for solving with a fixed time step from OdeBaseRK. It implements functions using an adaptive time step with different output options. Output values include only points along the solution, no interpolation (no dense output). The class handles FSAL methods properly if the fsal variable is `true` AND the "last" k vector that is the same as the first is identified in the klast variable (these variables live in OdeBaseRK). This is necessary to handle step rejection after k values have been changed during a step. klast should be set in a derived class. This class must also be initialized with the order of the lower order solution used for error estimation and it assumes that it will be using solutions of order errord and errord+1 for error estimation. That assumption means the class will not properly compute error estimates if the solution in yhi is not a single order higher than the solution in ylo.

This class is structured so that a derived class can define the stepping mechanism by implementing the step() function. Then, a grandchild class defines the system of equations to be solved by implementing the ode_funk() function.
*/

#ifndef ODE_BASE_ARK_HPP
#define ODE_BASE_ARK_HPP

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "OdeBaseRK.hpp"
#include "ode_misc.hpp"

class OdeBaseARK : public OdeBaseRK {
public:

    //variables for adaptive time stepping
    double dtmax; //hard upper limit for time step size, default is infinity
    double facsafe, facmin, facmax; //fractions for time step selection
    double abstol, reltol; //local absolute and relative error targets
    unsigned long long nrej; //counter for rejected steps

    //constructor
    OdeBaseARK (long int neq_, int errord_);
    //destructor
    ~OdeBaseARK ();
    //solve to a particular time with adaptive time stepping
    //    tend - stopping time, time at end of solution
    //    dt - time step
    //    snaps - number of snapshots to output
    //    dirout - relative path to preexisting directory for output
    void solve_adaptive (double tend, double dt0); //no output
    void solve_adaptive (double tend, double dt0, const char *dirout_); //full output
    void solve_adaptive (double tend, double dt0, int snaps, const char *dirout_); //snapshot output


protected:

    //order of the LOWER order solution used for error estimation
    int errord;

    //arrays for high and low order solutions in error estimation
    double *sollo, *solhi;
    //array for previous step's solution (for step reversal and maybe dense output)
    double *solprev;

    //reverse a cancelled step
    void reverse_step ();

    //functions for calculating the error at a step and optimal new time step
    double error (double abstol, double reltol);
    double facopt (double err);

};

#endif
