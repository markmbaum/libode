/*
This class provides a large vector containing the slope values of all stages with pointers to each of the individual stages. 
*/

#ifndef ODE_IRK_H_
#define ODE_IRK_H_

class OdeIRK {
    public:
        OdeIRK (unsigned long neq, int nk);
        ~OdeIRK ();
    protected:
        int nk_;
        double *kall_;
        double **k_;
};

#endif
