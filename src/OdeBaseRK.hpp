/*
This is a base class for Runge Kutta methods, with management of FSAL methods (adaptive or not) and implementation of three different versions of a function to solve a system of ODEs to a given time with a fixed time step (solve_fixed), each version generating different output or no output.

This class is structured so that a derived class can define the stepping mechanism by implementing the step() function. Then, a further derived class defines the system of equations to be solved by implementing the ode_funk() function.
*/

#ifndef ODE_BASE_RK_HPP
#define ODE_BASE_RK_HPP

#include <vector>
#include <string>

#include "ode_misc.hpp"
#include "OdeBase.hpp"

class OdeBaseRK : public OdeBase {
public:

    //constructor
    OdeBaseRK (unsigned long long neq_);
    //destructor
    ~OdeBaseRK ();

    //solve to a particular time with a fixed time step
    //    tend - stopping time, time at end of solution
    //    dt - time step
    //    snaps - number of snapshots to output
    //    dirout - relative path to preexisting directory for output
    void solve_fixed (double tend, double dt_); //no output, fast
    void solve_fixed (double tend, double dt_, const char *dirout_); //full output
    void solve_fixed (double tend, double dt_, int snaps, const char *dirout_); //some output

protected:

    //temporary y values for moving through stages
    double *soltemp;
    //flag for fsal methods
    bool fsal;
    //the last k values, MUST BE SET IN DERIVED FSAL METHODS IF NEEDED
    double *klast;

    /*function to see whether the next normal time step (not short step for snapping) will hit the termination time (tend)*/
    bool solve_done (double dt_, double tend);

};

#endif
