/*
A nonlinear diffusion equation in 1D
*/

#include <cmath>
#include <cstdio>

template<class Integrator>
class NonlinearDiffusion : public Integrator {
    public:
        NonlinearDiffusion (int nx_, double L_) : Integrator (nx_), nx(nx_), L(L_) {
            h = L/double(nx);
            d = new double[nx-1];
            for (int i=0; i<nx-1; i++)
                d[i] = INFINITY;
            xc = new double[nx];
            for (int i=0; i<nx; i++)
                xc[i] = h/2.0 + i*h;
        }
        ~NonlinearDiffusion () {
            delete [] d;
            delete [] xc;
        }

        //number of cells
        int nx;
        //domain length
        double L;
        //grid spacing
        double h;
        //diffusion exponent
        double gam;
        //cell centers
        double *xc;
        //diffusion coefficients
        double *d;

        //computes the diffusion coefficient
        double f_d (double u) { return( u*u*u ); }

        void ode_fun (double *solin, double *fout) {

            int i;
            double *u = solin;
            double q;

            //zero out fluxes before accumulating
            for (i=0; i<nx; i++) fout[i] = 0.0;

            //evaluate fluxes
            for (i=0; i<nx-1; i++) {
                //diffusion coefficient
                d[i] = f_d( (u[i] + u[i+1])/2.0 );
                //diffusive flux
                q = -d[i]*(u[i+1] - u[i])/h;
                //time derivatives
                fout[i]   -= q/h;
                fout[i+1] += q/h;
            }
        }

        double dt_adapt () {
            //compute maximum current diffusion equation
            double dmax = 0;
            for (int i=0; i<nx-1; i++)
                if ( d[i] > dmax )
                    dmax = d[i];
            //return the largest stable time step, with a little extra room
            return ( 0.45*(h*h/dmax) );
        }
};
