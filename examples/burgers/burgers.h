/*
Inviscid Burger's equation with a simple finite-volume upwinding scheme
*/

#include <cmath>
#define pi 3.1415926535897932384626433832795

#include <ode_linalg.h>

template<class Integrator>
class Burgers : public Integrator {
    public:
        Burgers (int nx_, double L_) : Integrator (nx_), nx(nx_), L(L_) {
            h = L/double(nx);
            xc = new double[nx];
            for (int i=0; i<nx; i++) xc[i] = h/2.0 + i*h;
        }
        ~Burgers () { delete [] xc; }

        //number of cells
        int nx;
        //domain length
        double L;
        //grid spacing
        double h;
        //cell centers
        double *xc;

        void ode_fun (double *solin, double *fout) {

            int i;
            double *u = solin;
            double f;

            //zero out fluxes before accumulating
            for (i=0; i<nx; i++) fout[i] = 0.0;

            //evaluate fluxes
            for (i=0; i<nx-1; i++) {
                //upwind
                f = (u[i+1] + u[i])/2.0 < 0.0 ? u[i+1]*u[i+1]/2.0 : u[i]*u[i]/2.0;
                fout[i]   -= f/h;
                fout[i+1] += f/h;
            }
        }

};
