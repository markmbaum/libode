#ifndef ODE_RK_H_
#define ODE_RK_H_

//!Provides space for stage slope values
class OdeRK {

    public:
        //!constructs
        /*!
        \param[in] neq size of ODE sytem
        \param[out] nk number of stages
        */
        OdeRK (unsigned long neq, int nk);
        //!destructs
        ~OdeRK ();

    protected:
        //!stage evaluations
        double **k_;
        //!number of stages, or k vectors
        int nk_;

};

#endif
