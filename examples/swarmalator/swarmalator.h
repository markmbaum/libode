#include <cstdlib>

#define pi 3.1415926535897932384626433832795

//------------------------------------------------------------------------------

//quick function for getting a pseudorandom number between 1 and 0
double rand_ () { return(double(rand())/double(RAND_MAX)); }

//------------------------------------------------------------------------------
/*
Oscillators that sync and swarm, or SWARMALATORS, from:

    3K. P. O’Keeffe, H. Hong, and S. H. Strogatz, Oscillators that sync and swarm, Nat. Commun. 8, 1504 (2017). doi:10.1038/s41467-017-01190-3

Evaluation of the system could be sped up by pre-calculating distances to avoid duplication, but oh well...
*/

template<class Integrator>
class Swarmalator : public Integrator {

    public:

        //constants
        double A, B, J, K;

        //each agent gets three equations, x, y, and phase
        Swarmalator (int n) : Integrator (3*n), n_(n), N_(double(n)) {

            //initial x coord, y coord, and phase (theta)
            double r1, r2;
            for (int i=0; i<n; i++) {
                //get random x,y coords in the unit circle
                r1 = 2*rand_() - 1;
                r2 = 2*rand_() - 1;
                while (sqrt(r1*r1 + r2*r2) > 1.0) {
                    r1 = 2*rand_() - 1;
                    r2 = 2*rand_() - 1;
                }
                Integrator::set_sol(i, r1);
                Integrator::set_sol(i+n, r2);
                //get a random phase in [0,2*pi)
                Integrator::set_sol(i+2*n, 2*pi*rand_());
            }
        }

        void ode_fun (double *solin, double *fout) {

            //alias
            double *x = solin,
                   *y = solin + n_,
                   *theta = solin + 2*n_;
            //variables for distances and such
            double delx, dely, cdist, tdist, dxdt, dydt, dtdt;

            //zero out the time derivatives before incrementing
            for (int i=0; i<3*n_; i++) fout[i] = 0.0;

            //inefficient calculation of time derivatives b/c distances are repeated
            for (int i=0; i<n_-1; i++) {
                for (int j=i+1; j<n_; j++) {
                    //distances
                    delx = x[j] - x[i];
                    dely = y[j] - y[i];
                    cdist = sqrt(delx*delx + dely*dely);
                    tdist = theta[j] - theta[i];
                    //time derivaties
                    dxdt = ((delx/cdist)*(A + J*cos(tdist))  - B*(delx/(cdist*cdist)))/N_;
                    fout[i]      += dxdt;
                    fout[j]      -= dxdt;
                    dydt = ((dely/cdist)*(A + J*cos(tdist))  - B*(dely/(cdist*cdist)))/N_;
                    fout[i+n_]   += dydt;
                    fout[j+n_]   -= dydt;
                    dtdt = (K*sin(tdist)/cdist)/N_;
                    fout[i+2*n_] += dtdt;
                    fout[j+2*n_] -= dtdt;
                }
            }
        }

    private:
        //number of agents
        int n_;
        double N_;
};
