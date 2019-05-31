/*
OdeBase is an abstract base class defining basic variables and functions used by methods of solving systems of ordinary differential equations (ODEs). It contains an array (sol) for the solution to each ODE at a given time, the time of the solution (t), and counters for function evaluations (neval) and steps (nstep). The ode_funk function defines the system of equations to solve. The step function should advance the solution forward a single time step, increment the time, increment neval, and increment nstep.
*/

#ifndef ODE_BASE_HPP
#define ODE_BASE_HPP

#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "ode_misc.hpp"

class OdeBase {
public:

    //number of equations in the system of ODEs
    unsigned long long neq;
    //time, initialized to zero
    double t;
    //time step
    double dt;
    //array for the solution, changing over time
    double *sol;
    //number of time steps
    unsigned long long nstep;
    //function evaluation counter, must be incremented in step() when defined
    unsigned long long neval;
    //interval of steps after which to check for nans and infs (zero to ignore)
    unsigned long long ncheck;

    //constructor
    OdeBase (unsigned long long neq_);
    //destructor
    ~OdeBase ();

    //the system of ODEs, defined in derived class
    virtual void ode_funk (double t, double *solin, double *fout) = 0;
    //advance a single time step, defined in derived class
    virtual void step (double dt) = 0;

    //write y values to file and handle snapshot variables as solution progresses
    void snap (std::string dirout, long *isnap, double *tsnap, double tend,
                int snaps, std::vector<double> *tout);

    //do any extra stuff before starting a solve (to be overridden in derived class)
    virtual void before_solve(std::string dirout);
    //do any extra stuff after stepping (to be overridden in derived class)
    //can only be used with solve_fixed() because adaptive solving sometimes
    //reverses steps
    virtual void step_extra (double t_);
    //do any extra stuff upon snapping (to be overridden in derived class)
    virtual void snap_extra (std::string dirout, long isnap, double tsnap);
    //do any extra stuff after completing a solve (to be overridden in derived class)
    virtual void after_solve(std::string dirout);

    //check solution for nans and infs
    void check_sol_integrity ();
};

#endif
