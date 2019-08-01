#ifndef ODE_RK_H_
#define ODE_RK_H_

class OdeRK {

    public:
        //constructor
        OdeRK (unsigned long neq, int nk);
        //destructor
        ~OdeRK ();

    protected:
        //stage evaluations
        double **k_;
        //number of stages, or k vectors
        int nk_;

};

#endif
