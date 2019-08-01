template<class Integrator>
class Lorentz : public Integrator {

    public:

        //parameters
        double sigma, beta, rho;

        //constructor
        Lorentz () : Integrator (3) {}

        //system of equations
        void ode_fun (double *solin, double *fout) {
            //alias
            double x = solin[0],
                   y = solin[1],
                   z = solin[2];
            //evaluate derivatives
            fout[0] = sigma*(y - x);          // dx/dt
            fout[1] = x*(rho - z) - y; // dy/dt
            fout[2] = x*y - beta*z ;   // dz/dt
        }
        //jacobian
        void ode_jac (double *solin, double **Jout) {
            //alias
            double x = solin[0],
                   y = solin[1],
                   z = solin[2];
            Jout[0][0] =  -sigma; Jout[0][1] = sigma; Jout[0][2] =   0.0;
            Jout[1][0] = rho - z; Jout[1][1] =  -1.0; Jout[1][2] =    -x;
            Jout[2][0] =       y; Jout[2][1] =     x; Jout[2][2] = -beta;
        }
};
