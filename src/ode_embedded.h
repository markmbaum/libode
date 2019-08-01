/*
This class implements the error analysis and time step step selection for embedded RK methods
*/

#ifndef ODE_EMBEDDED_H_
#define ODE_EMBEDDED_H_

#include <cmath>

#include "ode_util.h"
#include "ode_adaptive.h"

class OdeEmbedded : public OdeAdaptive {

    public:
        //constructor
        OdeEmbedded (unsigned long neq, bool need_jac, int lowerrord);
        //destructor
        ~OdeEmbedded ();

        //-------------------
        //getters and setters

        double get_facsafe () { return(facsafe_); }
        double get_facmin () { return(facmin_); }
        double get_facmax () { return(facmax_); }

        void set_facsafe (double facsafe) { facsafe_ = facsafe; }
        void set_facmin (double facmin) { facmin_ = facmin; }
        void set_facmax (double facmax) { facmax_ = facmax; }

    protected:
        //fractions for time step selection
        double facsafe_, facmin_, facmax_;

        //array for the lower order solution for error estimation
        //the higher order solution is assumed to go straight into sol
        double *solemb_; //embedded solution
        //calculates error estimate with lower and higher order solutions
        double error (double abstol, double reltol);
        //calculates factor for "optimal" next time step
        double facopt (double err);

        //------------------------------------------------
        //implementation of inherited methods for adapting

        //does the calculations to determine isrej and dtopt
        void adapt (double abstol, double reltol);
        //simply returns isrej
        bool is_rejected () { return( isrej_ ); }
        //simply returns dtopt
        double dt_next () { return( dtopt_ ); }

    private:
        //order of the LOWER order solution used for error estimation
        int lowerrord_;
        //flag for rejecting a step
        bool isrej_;
        //next time step
        double dtopt_;
};

#endif
