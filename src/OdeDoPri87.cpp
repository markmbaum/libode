#include "OdeDoPri87.hpp"

OdeDoPri87::OdeDoPri87 (unsigned long long neq_) : OdeBaseARK (neq_, 7) {

    //time stepping arrays
    k1 = new double[neq];
    k2 = new double[neq];
    k3 = new double[neq];
    k4 = new double[neq];
    k5 = new double[neq];
    k6 = new double[neq];
    k7 = new double[neq];
    k8 = new double[neq];
    k9 = new double[neq];
    k10 = new double[neq];
    k11 = new double[neq];
    k12 = new double[neq];
    k13 = new double[neq];
    //tableau of coefficents
    c2 =                   1.0/18; a21 =                   1.0/18;
    c3 =                   1.0/12; a31 =                   1.0/48; a32 = 1.0/16;
    c4 =                    1.0/8; a41 =                   1.0/32; a42 =    0.0; a43 =   3.0/32;
    c5 =                   5.0/16; a51 =                   5.0/16; a52 =    0.0; a53 = -75.0/64; a54 =                     75.0/64;
    c6 =                    3.0/8; a61 =                   3.0/80; a62 =    0.0; a63 =      0.0; a64 =                      3.0/16; a65 =                   3.0/20;
    c7 =                 59.0/400; a71 =     29443841.0/614563906; a72 =    0.0; a73 =      0.0; a74 =        77736538.0/692538347; a75 =   -28693883.0/1125000000; a76 =      23124283.0/1800000000;
    c8 =                 93.0/200; a81 =     16016141.0/946692911; a82 =    0.0; a83 =      0.0; a84 =        61564180.0/158732637; a85 =     22789713.0/633445777; a86 =     545815736.0/2771057229; a87 =    -180193667.0/1043307555;
    c9 =  5490023248.0/9719169821; a91 =     39632708.0/573591083; a92 =    0.0; a93 =      0.0; a94 =      -433636366.0/683701615; a95 =  -421739975.0/2616292301; a96 =      100302831.0/723423059; a97 =      790204164.0/839813087; a98 =     800635310.0/3783071287;
    c10 =                 13.0/20; a101 =  246121993.0/1340847787; a102 =   0.0; a103 =     0.0; a104 = -37695042795.0/15268766246; a105 = -309121744.0/1061227803; a106 =     -12992083.0/490766935; a107 =   6005943493.0/2108947869; a108 =    393006217.0/1396673457; a109 =    123872331.0/1001029789;
    c11 = 1201146811.0/1299019798; a111 = -1028468189.0/846180014; a112 =   0.0; a113 =     0.0; a114 =     8478235783.0/508512852; a115 = 1311729495.0/1432422823; a116 = -10304129995.0/1701304382; a117 = -48777925059.0/3047939560; a118 =  15336726248.0/1032824649; a119 = -45442868181.0/3398467696; a1110 =  3065993473.0/597172653;
    c12 =                     1.0; a121 =   185892177.0/718116043; a122 =   0.0; a123 =     0.0; a124 =    -3185094517.0/667107341; a125 = -477755414.0/1098053517; a126 =    -703635378.0/230739211; a127 =   5731566787.0/1027545527; a128 =    5232866602.0/850066563; a129 =   -4093664535.0/808688257; a1210 = 3962137247.0/1805957418; a1211 =   65686358.0/487910083;
    c13 =                     1.0; a131 =   403863854.0/491063109; a132 =   0.0; a133 =     0.0; a134 =    -5068492393.0/434740067; a135 =  -411421997.0/543043805; a136 =     652783627.0/914296604; a137 =   11173962825.0/925320556; a138 = -13158990841.0/6184727034; a139 =   3936647629.0/1978049680; a1310 =  -160528059.0/685178525; a1311 = 248638103.0/1413531060; a1312 =                   0.0;
                                   b1 =      14005451.0/335480064; b2 =     0.0; b3 =       0.0; b4 =                          0.0; b5 =                       0.0; b6 =      -59238493.0/1068277825; b7 =       181606767.0/758867731; b8 =       561292985.0/797845732; b9 =    -1041891430.0/1371343529; b10 =    760417239.0/1151165299; b11 =    118820643.0/751138087; b12 = -528747749.0/2220607170; b13 = 1.0/4;
                                   d1 =      13451932.0/455176632; d2 =     0.0; d3 =       0.0; d4 =                          0.0; d5 =                       0.0; d6 =      -808719846.0/976000145; d7 =     1757004468.0/5645159321; d8 =       656045339.0/265891186; d9 =    -3867574721.0/1518517206; d10 =     465885868.0/322736535; d11 =     53011238.0/667516719; d12 =                  2.0/45; d13   = 0.0;
}

OdeDoPri87::~OdeDoPri87 () {
    delete [] k1;
    delete [] k2;
    delete [] k3;
    delete [] k4;
    delete [] k5;
    delete [] k6;
    delete [] k7;
    delete [] k8;
    delete [] k9;
    delete [] k10;
    delete [] k11;
    delete [] k12;
    delete [] k13;
}

void OdeDoPri87::step (double dt) {

    unsigned long long i;

    //first rk step
    for (i=0; i<neq; i++) ode_funk(t, sol, k1);

    //k2
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*a21*k1[i];
    ode_funk(t + dt*c2, soltemp, k2);

    //k3
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a31*k1[i]
                                                  + a32*k2[i]);
    ode_funk(t + dt*c3, soltemp, k3);

    //k4
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a41*k1[i]
                                                  + a43*k3[i]);
    ode_funk(t + dt*c4, soltemp, k4);

    //k5
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a51*k1[i]
                                                  + a53*k3[i]
                                                  + a54*k4[i]);
    ode_funk(t + dt*c5, soltemp, k5);

    //k6
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a61*k1[i]
                                                  + a64*k4[i]
                                                  + a65*k5[i]);
    ode_funk(t + dt*c6, soltemp, k6);

    //k7
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a71*k1[i]
                                                  + a74*k4[i]
                                                  + a75*k5[i]
                                                  + a76*k6[i]);
    ode_funk(t + dt*c7, soltemp, k7);

    //k8
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a81*k1[i]
                                                  + a84*k4[i]
                                                  + a85*k5[i]
                                                  + a86*k6[i]
                                                  + a87*k7[i]);
    ode_funk(t + dt*c8, soltemp, k8);

    //k9
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a91*k1[i]
                                                  + a94*k4[i]
                                                  + a95*k5[i]
                                                  + a96*k6[i]
                                                  + a97*k7[i]
                                                  + a98*k8[i]);
    ode_funk(t + dt*c9, soltemp, k9);

    //k10
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a101*k1[i]
                                                  + a104*k4[i]
                                                  + a105*k5[i]
                                                  + a106*k6[i]
                                                  + a107*k7[i]
                                                  + a108*k8[i]
                                                  + a109*k9[i]);
    ode_funk(t + dt*c10, soltemp, k10);

    //k11
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a111*k1[i]
                                                  + a114*k4[i]
                                                  + a115*k5[i]
                                                  + a116*k6[i]
                                                  + a117*k7[i]
                                                  + a118*k8[i]
                                                  + a119*k9[i]
                                                  + a1110*k10[i]);
    ode_funk(t + dt*c11, soltemp, k11);

    //k12
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a121*k1[i]
                                                  + a124*k4[i]
                                                  + a125*k5[i]
                                                  + a126*k6[i]
                                                  + a127*k7[i]
                                                  + a128*k8[i]
                                                  + a129*k9[i]
                                                  + a1210*k10[i]
                                                  + a1211*k11[i]);
    ode_funk(t + dt*c12, soltemp, k12);

    //k13
    for (i=0; i<neq; i++) soltemp[i] = sol[i] + dt*(a131*k1[i]
                                                  + a134*k4[i]
                                                  + a135*k5[i]
                                                  + a136*k6[i]
                                                  + a137*k7[i]
                                                  + a138*k8[i]
                                                  + a139*k9[i]
                                                  + a1310*k10[i]
                                                  + a1311*k11[i]);
    ode_funk(t + dt*c13, soltemp, k13);

    //compute values for the whole step
    for (i=0; i<neq; i++) {
        //eighth order result
        solhi[i] = sol[i] + dt*(b1*k1[i]
                              + b6*k6[i]
                              + b7*k7[i]
                              + b8*k8[i]
                              + b9*k9[i]
                              + b10*k10[i]
                              + b11*k11[i]
                              + b12*k12[i]
                              + b13*k13[i]);
        //seventh order result
        sollo[i] = sol[i] + dt*(d1*k1[i]
                              + d6*k6[i]
                              + d7*k7[i]
                              + d8*k8[i]
                              + d9*k9[i]
                              + d10*k10[i]
                              + d11*k11[i]
                              + d12*k12[i]);
    }

    //take the 8th order solution
    for (i=0; i<neq; i++) {
        solprev[i] = sol[i];
        sol[i] = solhi[i];
    }

    //increment the time
    t += dt;
    //update function evaluation counter
    neval += 13;
    //add a step to the counter
    nstep++;
}
