# libode

This repo contains a collection of C++ and Python classes for solving systems of ordinary differential equations (ODEs) and a few programs to test the behavior of the different solvers. So far, only single-step, explicit methods are included in the C++ classes, but the Python classes include some low order implicit solvers. All of the solvers can be used with a fixed time step and several of them are embedded Runge-Kutta methods with adaptive time step selection. Some of them also take advantage of the first-same-as-last (FSAL) property. The C++ classes were originally modeled after Chris Rycroft's [example classes](https://github.com/chr1shr/am225_examples/tree/master/1a_ode_solvers) and have mostly been used to handle the temporal discretization in systems of PDEs (the [method of lines](https://en.wikipedia.org/wiki/Method_of_lines) approach). Except for examining the behavior of high-order solvers and providing a little extra flexibility, the Python classes are not more useful than scipy's [classes](https://docs.scipy.org/doc/scipy/reference/integrate.html). They have been useful, however, for testing some aspects of the solvers before implementing them in C++.

## Using the Solvers

#### C++

To use a C++ solver, a new C++ class must be created to inherit from one of the solver classes. This new inheriting class must
1. define the system of ODEs to be solved by implementing the `ode_funk` function. This is a virtual function in the base classes. Once it is implemented, it can be used by the stepping and solving functions in the base classes.
2. set initial conditions by setting values in the `sol` variable, which is an array with the same length as the system of ODEs.

Other than defining the system of equations and setting initial conditions, the derived class could also store whatever information and implement whatever other methods are necessary. This could be something simple like some extra functions for setting initial conditions, but it could be any other system that needs to run on top of an ODE solver, like the spatial discretization of a big PDE solver.

Some examples of simple systems of 2 ODEs, which have been used to test convergence and accuracy, are implemented [here](test/cpp/ode_explicit_test_systems.hpp) file. A template class is also provided [here](TemplateSolver.cpp) with a `main` function to drive it.

Each solver class has the `solve_fixed` function and, if it's an adaptive class, the `solve_adaptive` function. These functions return nothing and both have three call signatures.
* `void solve_fixed (double tend, double dt)`
* `void solve_fixed (double tend, double dt, const char *dirout)`
* `void solve_fixed (double tend, double dt, int snaps, const char *dirout)`

The first is a simple, bare version that advances the solution to a given point in time without any output. The system is integrated to `tend` with a time step of `dt` (or an initial step size `dt` for adaptive solving). The second call signature will write the complete solution for each variable into a specified directory as a binary file (full output). The third signature will write some number of evenly spaced snapshots into a specified directory. The signatures for `solve_adaptive` are identical, but `dt` is only an initial guess at the time step for adaptive solves.

The Makefile compiles all of the necessary code into the obj folder, then archives it in the bin directory as a file called libode.a. Test programs are compiled with `make tests` and they can all be run at once with the `run_tests.sh` script. Before compiling, the `_config.mk` file must be copied and renamed `config.mk`. In that file, the desired compiler and other settings are indicated. Finally, to use the solvers, you could point a compiler toward `libode.a` or toward the object files directly, in addition to the header files in the `src` folder.

#### Python

The Python classes are structured much like the C++ ones and they have the same `solve_fixed` and `solve_adaptive` methods. To access them, you need to have the main repository folder in your path, then use `from pylibode import solvers` to import a DataFrame of all of the individual solvers. The DataFrame's rows are indexed by the short name of each solver (Euler, RK4, ...) and has each class's constructor in one of its columns. Each constructor requires a function defining the system of ODEs and an iterable of the initial values. For example, to use Euler's method from the DataFrame, you would do something like
```
from pylibode import solvers
s = solvers['constructor']['Euler'](ode_funk, [1., 2., 3.])
t, y = s.solve_fixed(10, 1e-3)
```
which would solve the system of ODEs defined by `ode_funk`, with initial conditions (1,2,3), for 10 units of the independent variable, and a step size of 0.001. Alternatively, you could import the solver class directly with `from pylibode import OdeEuler`, then use it without the DataFrame.

## Testing

The convergence and accuracy of all the C++ solvers has been tested using a few programs with source files named `test/ode_test_*` and executables named `bin/ode_test_*`. Python scripts to plot the output of these programs are named `scripts/plot_*`. All these tests can be compiled, run, and plotted with the `run_tests.sh` script. Python solvers can be tested with scripts in the `test/python` directory.

The image below shows the results of the adaptive fifth order method DOPRI5 (Dormand-Prince 5), which is the method used by MATLAB's `ode45` function, solving a system of two coupled oscillators. The plot was generated from the top directory with:
```
> make tests
> ./bin/test_adapt.exe
> cd scripts
> python plot_adapt.py
```
![-](img/dopri54_fsal_osc2_adapt.png)

The `test_work` program examines how many function evaluations (the number of times the system of ODEs is evaluated) are required for different levels of accuracy for several methods. Each method solves an identical initial value problem and the final values are compared to a known analytical solution. The plot below, generated with `scripts/plot_work.py`, shows the results. As expected, but still amazingly, the eighth order method of Dormand and Prince (`OdeDoPri87`) achieves error reduction by a factor of about 100,000,000 when the number of function calls is increased by a factor of 10.

![-](img/work.png)

Another test program, `test_conv`, can be used to confirm the order of accuracy of different methods. Below are the results for a strong stability preserving method of order three. The black dots fall neatly along a line proportional to the cube of the step size until machine precision is reached for the global error, confirming the third order accuracy.

![-](img/ssp3_conv.png)
