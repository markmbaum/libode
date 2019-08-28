#ifndef ODE_EMBEDDED_H_
#define ODE_EMBEDDED_H_

//! \file ode_embedded.h

#include <cmath>

#include "ode_util.h"
#include "ode_adaptive.h"

//!Base clase implementing methods for embedded Runge-Kutta error estimation
class OdeEmbedded : public OdeAdaptive {

    public:
        //!constructs
        /*!
        \param[in] neq size of ODE system
        \param[in] need_jac flag signaling whether the Jacobian of the system is needed
        \param[in] lowerord the order of the lower order embedded method (assumed to be one less than the higher order one)
        */
        OdeEmbedded (unsigned long neq, bool need_jac, int lowerord);
        //!destructs
        ~OdeEmbedded ();

        //-------------------
        //getters and setters

        //!gets safety factor applied to time step selection
        double get_facsafe () { return(facsafe_); }
        //!gets minimum allowable fraction change in time step
        double get_facmin () { return(facmin_); }
        //!gets maximum allowable fraction change in time step
        double get_facmax () { return(facmax_); }

        //!sets safety factor applied to time step selection
        void set_facsafe (double facsafe) { facsafe_ = facsafe; }
        //!sets minimum allowable fraction change in time step
        void set_facmin (double facmin) { facmin_ = facmin; }
        //!sets maximum allowable fraction change in time step
        void set_facmax (double facmax) { facmax_ = facmax; }

    protected:

        //!safety factor applied to time step selection
        double facsafe_;
        //!minimum allowable fraction change in time step
        double facmin_;
        //!maximum allowable fraction change in time step
        double facmax_;

        //!embedded solution array
        double *solemb_;
        //!calculates error estimate with lower and higher order solutions
        double error (double abstol, double reltol);
        //!calculates factor for "optimal" next time step
        double facopt (double err);

        //------------------------------------------------
        //implementation of inherited methods for adapting

        //!does the calculations to determine isrej and dtopt
        void adapt (double abstol, double reltol);
        //!simply returns isrej
        bool is_rejected () { return( isrej_ ); }
        //!simply returns dtopt
        double dt_next () { return( dtopt_ ); }

    private:
        //!order of the LOWER order solution used for error estimation
        int lowerord_;
        //!flag for rejecting a step
        bool isrej_;
        //!next time step
        double dtopt_;
};

#endif
