#include <cmath>
#include <vector>

#include "ode_io.h"

#define pi 3.14159265358979323846264338327950288419716939937510582


namespace ode {

template<class Integrator>
class DoublePendulum : public Integrator {
    public:
        //physical parameters
        double g,  //gravity (m/s^2)
               m1, //mass 1 (kg)
               m2, //mass 2 (kg)
               L1, //length 1 (m)
               L2; //length 2 (m)

        //position vectors
        std::vector<double> x1, y1, x2, y2;

        //constructor
        DoublePendulum () : Integrator (4) {
            for (int i=0; i<4; i++) Integrator::set_sol(i, 0.0);
        }

        //system of equations
        void ode_fun (double *solin, double *fout) {

            //http://scienceworld.wolfram.com/physics/DoublePendulum.html

            //alias
            double t1 = solin[0],
                   p1 = solin[1],
                   t2 = solin[2],
                   p2 = solin[3];
            //evaluate
            double c1, c2;
            c1 = p1*p2*sin(t1 - t2)/(L1*L2*(m1 + m2*sin(t1 - t2)*sin(t1 - t2)));
            c2 = (L2*L2*m2*p1*p1 + L1*L1*(m1 + m2)*p2*p2 - L1*L2*m2*p1*p2*cos(t1 - t2))*sin(2*(t1 - t2))/(2*L1*L1*L2*L2*(m1 + m2*sin(t1 - t2)*sin(t1 - t2)));
            fout[0] = (L2*p1 - L1*p2*cos(t1 - t2))/(L1*L1*L2*(m1 + m2*sin(t1 - t2)*sin(t1 - t2)));
            fout[1] = -(m1 + m2)*g*L1*sin(t1) - c1 + c2;
            fout[2] = (L1*(m1 + m2)*p2 - L2*m2*p1*cos(t1 - t2))/(L1*L2*L2*m2*(m1 + m2*sin(t1 - t2)*sin(t1 - t2)));
            fout[3] = -m2*g*L2*sin(t2) + c1 - c2;
        }

    private:

        double f_x1 () {
            return(L1*sin(Integrator::get_sol(0)));
        }
        double f_y1 () {
            return(-L1*cos(Integrator::get_sol(0)));
        }
        double f_x2 () {
            return(f_x1() + L2*sin(Integrator::get_sol(2)));
        }
        double f_y2 () {
            return(f_y1() - L2*cos(Integrator::get_sol(2)));
        }

        void before_solve () {
            x1.clear(); x1.push_back( f_x1() );
            y1.clear(); y1.push_back( f_y1() );
            x2.clear(); x2.push_back( f_x2() );
            y2.clear(); y2.push_back( f_y2() );
        }
        void after_capture (double t) {
            (void)t;
            x1.push_back( f_x1() );
            y1.push_back( f_y1() );
            x2.push_back( f_x2() );
            y2.push_back( f_y2() );
        }
        void after_solve () {
            ode_write("out/x1", x1.data(), x1.size());
            ode_write("out/y1", y1.data(), y1.size());
            ode_write("out/x2", x2.data(), x2.size());
            ode_write("out/y2", y2.data(), y2.size());
        }
};

} // namespace ode
