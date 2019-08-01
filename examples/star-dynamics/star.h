#include <cmath>
#include <vector>

#include "ode_io.h"

template<class Integrator>
class Star : public Integrator {
    public:
        //non-dimensionalized physical parameters
        double a, b, c, A, C, Omega;
        //alias the solution array with position and momentum pointers
        double *q;
        double *p;
        //for storing the Hamiltonian along the solution
        std::vector<double> H;

        //constructor
        Star () : Integrator (6) {
            q = Integrator::get_sol();
            p = Integrator::get_sol() + 3;
        }

        //Hamiltonian
        double Hamiltonian (double *q, double *p) {
            return( (p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/2.0 + Omega*(p[0]*q[1] - p[1]*q[0]) + V(q) );
        }

        //system of equations
        void ode_fun (double *solin, double *fout) {
            double *q = solin,
                   *p = solin + 3;
            double dV = A/(C + q[0]*q[0]/(a*a) + q[1]*q[1]/(b*b) + q[2]*q[2]/(c*c));
            fout[0] = p[0] + Omega*q[1];
            fout[1] = p[1] - Omega*q[0];
            fout[2] = p[2];
            fout[3] = Omega*p[1] - dV*2*q[0]/(a*a);
            fout[4] = -Omega*p[0] - dV*2*q[1]/(b*b);
            fout[5] = -dV*2*q[2]/(c*c);
        }

    private:

        double V (double *q) {
            return( A*log(C + q[0]*q[0]/(a*a) + q[1]*q[1]/(b*b) + q[2]*q[2]/(c*c)) );
        }

        void before_solve () { H.clear(); H.push_back( Hamiltonian(q, p) ); }
        void after_step (double t) { (void)t; H.push_back( Hamiltonian(q, p) ); }
};
