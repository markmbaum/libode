/*
This class simply provides a temporary solution array for moving through RK stages with explicit solvers
*/

#ifndef ODE_ERK_H_
#define ODE_ERK_H_

class OdeERK {

    public:
        //constructor
        OdeERK (unsigned long neq);
        //destructor
        ~OdeERK ();

    protected:
        //temporary solution vector
        double *soltemp_;
};

#endif
