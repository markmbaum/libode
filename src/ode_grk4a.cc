//! \file ode_grk4a.cc

#include "ode_grk4a.h"

OdeGRK4A::OdeGRK4A (unsigned long neq) :
    OdeEmbedded (neq, true, 3),
    OdeRosenbrock (neq, 4) {

    method_ = "GRK4A";
    //main diagonal, gam_ii
    gam = 0.395;
    //all other gamma values
    gam21 = -0.767672395484;
    gam31 = -0.851675323742; gam32 =  0.522967289188;
    gam41 =  0.288463109545; gam42 = 0.0880214273381; gam43 = -0.337389840627;
    //reduced gamma values
    g21 = gam21/gam;
    g31 = gam31/gam; g32 = gam32/gam;
    g41 = gam41/gam; g42 = gam42/gam; g43 = gam43/gam;
    //alpha values
    alp21 = 0.438;
    alp31 = 0.796920457938; alp32 = 0.0730795420615;
    //k weights for 4th and 3rd order solutions
    b1 = 0.199293275701; b2 = 0.482645235674; b3 = 0.0680614886256; b4 = 0.25;
    d1 = 0.346325833758; d2 = 0.285693175712; d3 =  0.367980990530;

}

void OdeGRK4A::step_ (double dt) {

    unsigned long i;

    //------------------------------------------------------------------
    //first evaluate the Jacobian, do arithmetic, and LU factorize it

    ode_jac_(sol_, Jac_);
    prep_jac(Jac_, neq_, dt, p_);

    //------------------------------------------------------------------
    //stage 1

    ode_fun_(sol_, k_[0]);
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[0][i];
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[0]);

    //------------------------------------------------------------------
    //stage 2

    //temporary solution vector for evaluating ode_fun
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + alp21*k_[0][i];
    //put new slope values into k vector temporarily
    ode_fun_(soltemp_, k_[1]);
    //compute right hand side of matrix equation
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[1][i] + g21*k_[0][i];
    //solve the matrix equation
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[1]);
    //recover the actual k values
    for (i=0; i<neq_; i++) k_[1][i] -= g21*k_[0][i];

    //------------------------------------------------------------------
    //stage 3

    //temporary solution vector for evaluating ode_fun
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + alp31*k_[0][i]
                                                 + alp32*k_[1][i];
    //put new slope values into k vector temporarily
    ode_fun_(soltemp_, k_[2]);
    //compute right hand side of matrix equation
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[2][i] + g31*k_[0][i]
                                                 + g32*k_[1][i];
    //solve the matrix equation
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[2]);
    //recover the actual k values
    for (i=0; i<neq_; i++) k_[2][i] -= g31*k_[0][i] + g32*k_[1][i];

    //------------------------------------------------------------------
    //stage 4

    //temporary solution vector for evaluating ode_fun
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + alp31*k_[0][i]
                                                 + alp32*k_[1][i];
    //put new slope values into k vector temporarily
    ode_fun_(soltemp_, k_[3]);
    //compute right hand side of matrix equation
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[3][i] + g41*k_[0][i]
                                                 + g42*k_[1][i]
                                                 + g43*k_[2][i];
    //solve the matrix equation
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[3]);
    //recover the actual k values
    for (i=0; i<neq_; i++) k_[3][i] -= g41*k_[0][i] + g42*k_[1][i] + g43*k_[2][i];

    //------------------------------------------------------------------
    //solution and embedded solution

    for (i=0; i<neq_; i++) {
        solemb_[i] = sol_[i] + (d1*k_[0][i]
                              + d2*k_[1][i]
                              + d3*k_[2][i]);
        sol_[i]    = sol_[i] + (b1*k_[0][i]
                              + b2*k_[1][i]
                              + b3*k_[2][i]
                              + b4*k_[3][i]);
    }
}
