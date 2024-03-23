/*
Viscid Burger's equation with a simple finite-volume upwinding scheme
*/

#include <cmath>

namespace ode {
    
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
            double f,g;

            //zero out fluxes before accumulating
            for (i=0; i<nx; i++) fout[i] = 0.0;

            //evaluate fluxes
            for (i=0; i<nx-1; i++) {
                //upwind advective flux
                f = (u[i+1] + u[i])/2.0 < 0.0 ? u[i+1]*u[i+1]/2.0 : u[i]*u[i]/2.0;
                //diffusive flux
                g = -(u[i+1] - u[i])/(1000*h);
                //time derivatives
                fout[i]   -= (f + g)/h;
                fout[i+1] += (f + g)/h;
            }
        }
};

} // namespace ode
