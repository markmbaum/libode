#ifndef ODE_SYSTEMS_H_
#define ODE_SYSTEMS_H_

#include <cmath>

#include "ode_euler.h"
#include "ode_trapz.h"
#include "ode_ssp_3.h"
#include "ode_rkf_32.h"
#include "ode_rk_4.h"
#include "ode_rk_43.h"
#include "ode_rkck.h"
#include "ode_dopri_54.h"
#include "ode_vern_65.h"
#include "ode_vern_76.h"
#include "ode_dopri_87.h"
#include "ode_vern_98.h"
#include "ode_grk4a.h"
#include "ode_row6a.h"
#include "ode_backward_euler.h"
#include "ode_gauss_6.h"
#include "ode_lobatto_iiic_6.h"
#include "ode_radau_iia_5.h"
#include "ode_geng_5.h"
#include "ode_sdirk_43.h"

//Dahlquist test equation, the solution of which is just an exponential
template<class Integrator>
class Dahl : public Integrator {
    public:
        Dahl () : Integrator (2) {
            Integrator::set_name("Dahl");
            Integrator::set_sol(0, 1.0);
            Integrator::set_sol(1, 1.0);
        }
        void ode_fun (double *solin, double *fout) {
            fout[0] = -solin[0];
            fout[1] = -2*solin[1];
        }
        void ode_jac (double *solin, double **Jout) {
            (void)solin;
            Jout[0][0] = -1.0;
            Jout[0][1] =  0.0;
            Jout[1][0] =  0.0;
            Jout[1][1] = -2.0;
        }
};

//This is a system of two oscillating equations, which become increasingly
//steep with time. The solutions are
//  y1 = cos(t^2/2),   y1(0) = 1
//  y2 = sin(t^2/2),   y2(0) = 0
template<class Integrator>
class Osc1 : public Integrator {
    public:
        Osc1 () : Integrator (3) {
            Integrator::set_name("Osc1");
            Integrator::set_sol(0, 1.0);
            Integrator::set_sol(1, 0.0);
            Integrator::set_sol(2, 0.0); //time
        }
        void ode_fun (double *solin, double *fout) {
            double x = solin[0], y = solin[1], t = solin[2];
            fout[0] = -t*y;
            fout[1] = t*x;
            fout[2] = 1.0; //time
        }
        void ode_jac (double *solin, double **Jout) {
            double x = solin[0], y = solin[1], t = solin[2];
            Jout[0][0] = 0.0; Jout[0][1] =  -t; Jout[0][2] =  -y;
            Jout[1][0] =   t; Jout[1][1] = 0.0; Jout[1][2] =   x;
            Jout[2][0] = 0.0; Jout[2][1] = 0.0; Jout[2][2] = 0.0;
        }
};

//This is a second system of two oscillating equations, taken from Hairer II.4
//The solutions are:
//  y1 = exp(sin(t^2)),   y1(0) = 1
//  y2 = exp(cos(t^2)),   y2(0) = e
template<class Integrator>
class Osc2 : public Integrator {
    public:
        Osc2 () : Integrator (3) {
            Integrator::set_name("Osc2");
            Integrator::set_sol(0, 1.0);
            Integrator::set_sol(1, exp(1.0));
            Integrator::set_sol(2, 0.0); //time
        }
        void ode_fun (double *solin, double *fout) {
            solin[1] > 1e-3 ? fout[0] =  2*solin[2]*solin[0]*log(solin[1])
                            : fout[0] =  2*solin[2]*solin[0]*log(1e-3);
            solin[0] > 1e-3 ? fout[1] = -2*solin[2]*solin[1]*log(solin[0])
                            : fout[1] = -2*solin[2]*solin[1]*log(1e-3);
            fout[2] = 1.0; //time
        }
};

//This is a classic test system of two coupled, nonlinear equations called
//the Brusselator (Brussels oscillator).
//  https://en.wikipedia.org/wiki/Brusselator
template<class Integrator>
class Brus : public Integrator {
    public:
        Brus () : Integrator (2) {
            Integrator::set_name("Brus");
            Integrator::set_sol(0, 1.0);
            Integrator::set_sol(1, 1.0);
        }
        void ode_fun (double *solin, double *fout) {
            double x = solin[0], y = solin[1];
            fout[0] = 1.0 + x*x*y - 4.0*x;
            fout[1] = 3.0*x - x*x*y;
        }
        void ode_jac (double *solin, double **Jout) {
            double x = solin[0], y = solin[1];
            Jout[0][0] = 2*x*y - 4.0; Jout[0][1] = x*x;
            Jout[1][0] = 3.0 - 2*x*y; Jout[1][1] = -x*x;
        }
};

//The Van der Pol oscillator, a classic stiff problem
template<class Integrator>
class Vdp : public Integrator {
    public:
        //damping parameter
        double mu;
        Vdp () : Integrator(2) {
            mu = 100.0;
            Integrator::set_name("Vdp");
            Integrator::set_sol(0, 2.1);
            Integrator::set_sol(1, 0.0);
        }
        void ode_fun (double *solin, double *fout) {
            //alias
            double y1 = solin[0], y2 = solin[1];
            //compute derivatives
            fout[0] = y2;
            fout[1] = mu*mu*((1.0 - y1*y1)*y2 - y1);
        }
        void ode_jac (double *solin, double **Jout) {
            //alias
            double y1 = solin[0], y2 = solin[1];
            Jout[0][0] = 0.0;
            Jout[0][1] = 1.0;
            Jout[1][0] = -2.0*mu*mu*y1*y2 - mu*mu;
            Jout[1][1] = mu*mu*(1 - y1*y1);
        }
};

//stiff system from Hairer and Wanner
template<class Integrator>
class StiffCliff : public Integrator {
    public:
        StiffCliff () : Integrator(3) {
            Integrator::set_name("StiffCliff");
            Integrator::set_sol(0, 1.0);
            Integrator::set_sol(1, 0.0);
            Integrator::set_sol(2, 0.0);
            fac = 1000.0;
        }
        void ode_fun (double *solin, double *fout) {
            double y1 = solin[0], y2 = solin[1], x = solin[2];
            fout[0] = -fac*(cos(x)*y1 + sin(x)*y2 + 1.0);
            fout[1] = -fac*(-sin(x)*y1 + cos(x)*y2 + 1.0);
            fout[2] = 1.0;

        }
        void ode_jac (double *solin, double **Jout) {
            double y1 = solin[0], y2 = solin[1], x = solin[2];
            Jout[0][0] = -fac*cos(x); Jout[0][1] = -fac*sin(x); Jout[0][2] = -fac*(-sin(x)*y1 + cos(x)*y2);
            Jout[1][0] =  fac*sin(x); Jout[1][1] = -fac*cos(x); Jout[1][2] = -fac*(-cos(x)*y1 - sin(x)*y2);
            Jout[2][0] =         0.0; Jout[2][1] =         0.0; Jout[2][2] =                           0.0;
        }
    private:
        double fac;
};

#endif
