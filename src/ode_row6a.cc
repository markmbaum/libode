#include "ode_row6a.h"

OdeROW6A::OdeROW6A (unsigned long neq) :
    OdeBase (neq, true),
    OdeRosenbrock (neq, 6) {

    method_ = "ROW6A";
    //main diagonal, gam_ii
    gam = 0.33414236706805043;

    a21 =  0.66828473413610087;
    a31 =   0.5852480389573658; a32 = -0.048594008221492802;
    a41 = -0.61719233202999775; a42 =  -0.83995264476522158; a43 =   0.62641917900148600;
    a51 =   3.5406887484552165; a52 =   0.65991497772646308; a53 =  -0.63661180895697222; a54 = -1.1945984675295562;
    a61 =  0.80783664328582613; a62 =   0.10194631616818569; a63 = -0.078396778850607012; a64 = -0.044341977375427388; a65 = 0.013074732797453325;

    c21 = -5.8308828523185086;
    c31 = -4.0175939515896193; c32 =  0.43970131925236112;
    c41 =  7.7228006257490299; c42 =   4.3368108251435758; c43 = -2.8219574578033366;
    c51 = -1.0516225114542007; c52 = -0.58853585181331353; c53 =  2.0433794587212771; c54 = 5.0098631723809151;
    c61 = -6.7357785372199458; c62 = -0.53593889506199845; c63 = 0.38622517020810987; c64 = 0.21066472713931598; c65 = -0.053546655670373728;

    m1 = 11.358660043232931; m2 = -6.9896898855829058; m3 = -4.5967580421042947; m4 = -3.7220984696531517; m5 = 0.96012685868421520; m6 = 12.953396234292936;
}

void OdeROW6A::step_ (double dt) {

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
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + a21*k_[0][i];
    //put new slope values into k vector temporarily
    ode_fun_(soltemp_, k_[1]);
    //compute right hand side of matrix equation
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[1][i] + c21*k_[0][i];
    //solve the matrix equation
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[1]);

    //------------------------------------------------------------------
    //stage 3

    //temporary solution vector for evaluating ode_fun
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + a31*k_[0][i]
                                                 + a32*k_[1][i];
    //put new slope values into k vector temporarily
    ode_fun_(soltemp_, k_[2]);
    //compute right hand side of matrix equation
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[2][i] + c31*k_[0][i]
                                                 + c32*k_[1][i];
    //solve the matrix equation
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[2]);

    //------------------------------------------------------------------
    //stage 4

    //temporary solution vector for evaluating ode_fun
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + a41*k_[0][i]
                                                 + a42*k_[1][i]
                                                 + a43*k_[2][i];
    //put new slope values into k vector temporarily
    ode_fun_(soltemp_, k_[3]);
    //compute right hand side of matrix equation
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[3][i] + c41*k_[0][i]
                                                 + c42*k_[1][i]
                                                 + c43*k_[2][i];
    //solve the matrix equation
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[3]);

    //------------------------------------------------------------------
    //stage 5

    //temporary solution vector for evaluating ode_fun
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + a51*k_[0][i]
                                                 + a52*k_[1][i]
                                                 + a53*k_[2][i]
                                                 + a54*k_[3][i];
    //put new slope values into k vector temporarily
    ode_fun_(soltemp_, k_[4]);
    //compute right hand side of matrix equation
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[4][i] + c51*k_[0][i]
                                                 + c52*k_[1][i]
                                                 + c53*k_[2][i]
                                                 + c54*k_[3][i];
    //solve the matrix equation
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[4]);

    //------------------------------------------------------------------
    //stage 6

    //temporary solution vector for evaluating ode_fun
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + a61*k_[0][i]
                                                 + a62*k_[1][i]
                                                 + a63*k_[2][i]
                                                 + a64*k_[3][i]
                                                 + a65*k_[4][i];
    //put new slope values into k vector temporarily
    ode_fun_(soltemp_, k_[5]);
    //compute right hand side of matrix equation
    for (i=0; i<neq_; i++) rhs_[i] = dt*k_[5][i] + c61*k_[0][i]
                                                 + c62*k_[1][i]
                                                 + c63*k_[2][i]
                                                 + c64*k_[3][i]
                                                 + c65*k_[4][i];
    //solve the matrix equation
    ode_solve_LU(Jac_, p_, rhs_, neq_, k_[5]);

    //------------------------------------------------------------------
    //solution and embedded solution

    for (i=0; i<neq_; i++)
        sol_[i] = sol_[i] + (m1*k_[0][i]
                           + m2*k_[1][i]
                           + m3*k_[2][i]
                           + m4*k_[3][i]
                           + m5*k_[4][i]
                           + m6*k_[5][i]);
}
