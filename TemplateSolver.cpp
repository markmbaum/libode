#include <cstdio>
#include <cmath>

//#include "ode_misc.hpp"
//#include "OdeEuler.hpp"
//#include "OdeMidpnt.hpp"
//#include "OdeTrapz.hpp"
//#include "OdeSsp3.hpp"
//#include "OdeRKF23.hpp"
//#include "OdeRK4.hpp"
//#include "OdeRK34.hpp"
//#include "OdeRKCK.hpp"
//#include "OdeDoPri45.hpp"
//#include "OdeButcher6.hpp"
//#include "OdeDoPri78.hpp"

class TemplateSolver : public OdeSolver {
public:

    //solver doesn't have to take the number of equations as an input,
    //can be preset
    TemplateSolver (nequations) : OdeSolver (nequations) {
        //set initial conditions
    };

    ~TemplateSolver () {/* anything to do? */};

    //define the system of equations
    void ode_funk (double t_, double *solin, double *fout) {
        //solin contains the current state of the solution
        //fout is the vector of derivatives, populated by this function
    }
};

int main () {

    //instantiate the solver
    sys = TemplateSolver(/*nequations*/);

    //solve the system with different output options
    //sys.solve_fixed(tend, dt);
    //sys.solve_fixed(tend, dt, dirout);
    //sys.solve_fixed(tend, dt, nsnaps, dirout);

    printf("%lu function evaluations, %lu steps\n", sys.neval, sys.nstep);

    //solve the system with different output options
    //sys.solve_adaptive(tend, dt0);
    //sys.solve_adaptive(tend, dt0, dirout);
    //sys.solve_adaptive(tend, dt0, nsnaps, dirout);

    printf("%lu function evaluations, %lu steps, %lu rejected steps\n",
        sys.neval, sys.nstep, sys.nrej);

    return(0);
}
