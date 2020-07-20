//! \file ode_vern_76.cc

#include "ode_vern_76.h"

OdeVern76::OdeVern76 (unsigned long neq) :
    OdeEmbedded (neq, false, 6),
    OdeRK (neq, 10),
    OdeERK (neq) {

    method_ = "Vern76";

    /*
    coefficients copied from:
        http://people.math.sfu.ca/~jverner/
        http://people.math.sfu.ca/~jverner/RKV76.IIa.Efficient.00001675585.081206.CoeffsOnlyFLOAT
    */

    c2 =  0.005;
    c3 =  0.1088888888888888888888888888888888888889;
    c4 =  0.1633333333333333333333333333333333333333;
    c5 =  0.4555000000000000000000000000000000000000;
    c6 =  0.6095094489978381317087004421486024949638;
    c7 =  0.884;
    c8 =  0.925;
    c9 =  1.0;
    c10 = 1.0;

    a21 =  c2;

    a31 = -1.076790123456790123456790123456790123457;
    a32 =  1.185679012345679012345679012345679012346;

    a41 =  0.04083333333333333333333333333333333333333;
    a43 =  0.1225;

    a51 =  0.6389139236255726780508121615993336109954;
    a53 = -2.455672638223656809662640566430653894211;
    a54 =  2.272258714598084131611828404831320283215;

    a61 = -2.661577375018757131119259297861818119279;
    a63 =  10.80451388645613769565396655365532838482;
    a64 = -8.353914657396199411968048547819291691541;
    a65 =  0.8204875949566569791420417341743839209619;

    a71 =  6.067741434696770992718360183877276714679;
    a73 = -24.71127363591108579734203485290746001803;
    a74 =  20.42751793078889394045773111748346612697;
    a75 = -1.906157978816647150624096784352757010879;
    a76 =  1.006172249242068014790040335899474187268;

    a81 =  12.05467007625320299509109452892778311648;
    a83 = -49.75478495046898932807257615331444758322;
    a84 =  41.14288863860467663259698416710157354209;
    a85 = -4.461760149974004185641911603484815375051;
    a86 =  2.042334822239174959821717077708608543738;
    a87 = -0.09834843665406107379530801693870224403537;

    a91 =  10.13814652288180787641845141981689030769;
    a93 = -42.64113603171750214622846006736635730625;
    a94 =  35.76384003992257007135021178023160054034;
    a95 = -4.348022840392907653340370296908245943710;
    a96 =  2.009862268377035895441943593011827554771;
    a97 =  0.3487490460338272405953822853053145879140;
    a98 = -0.2714390051048312842371587140910297407572;

    a101 = -45.03007203429867712435322405073769635151;
    a103 =  187.3272437654588840752418206154201997384;
    a104 = -154.0288236935018690596728621034510402582;
    a105 =  18.56465306347536233859492332958439136765;
    a106 = -7.141809679295078854925420496823551192821;
    a107 =  1.308808578161378625114762706007696696508;

    b1 =  0.04715561848627222170431765108838175679569;
    b4 =  0.2575056429843415189596436101037687580986;
    b5 =  0.2621665397741262047713863095764527711129;
    b6 =  0.1521609265673855740323133199165117535523;
    b7 =  0.4939969170032484246907175893227876844296;
    b8 = -0.2943031171403250441557244744092703429139;
    b9 =  0.08131747232495109999734599440136761892478;

    d1 =  0.04460860660634117628731817597479197781432;
    d4 =  0.2671640378571372680509102260943837899738;
    d5 =  0.2201018300177293019979715776650753096323;
    d6 =  0.2188431703143156830983120833512893824578;
    d7 =  0.2289871705411202883378173889763552365362;
    d10 = 0.02029518466335628222767054793810430358554;
}

void OdeVern76::step_ (double dt) {

    unsigned long i;

    //------------------------------------------------------------------
    ode_fun_(sol_, k_[0]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*a21*k_[0][i];
    ode_fun_(soltemp_, k_[1]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a31*k_[0][i]
                                                     + a32*k_[1][i]);
    ode_fun_(soltemp_, k_[2]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a41*k_[0][i]
                                                     + a43*k_[2][i]);
    ode_fun_(soltemp_, k_[3]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a51*k_[0][i]
                                                     + a53*k_[2][i]
                                                     + a54*k_[3][i]);
    ode_fun_(soltemp_, k_[4]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a61*k_[0][i]
                                                     + a63*k_[2][i]
                                                     + a64*k_[3][i]
                                                     + a65*k_[4][i]);
    ode_fun_(soltemp_, k_[5]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a71*k_[0][i]
                                                     + a73*k_[2][i]
                                                     + a74*k_[3][i]
                                                     + a75*k_[4][i]
                                                     + a76*k_[5][i]);
    ode_fun_(soltemp_, k_[6]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a81*k_[0][i]
                                                     + a83*k_[2][i]
                                                     + a84*k_[3][i]
                                                     + a85*k_[4][i]
                                                     + a86*k_[5][i]
                                                     + a87*k_[6][i]);
    ode_fun_(soltemp_, k_[7]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a91*k_[0][i]
                                                     + a93*k_[2][i]
                                                     + a94*k_[3][i]
                                                     + a95*k_[4][i]
                                                     + a96*k_[5][i]
                                                     + a97*k_[6][i]
                                                     + a98*k_[7][i]);
    ode_fun_(soltemp_, k_[8]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) soltemp_[i] = sol_[i] + dt*(a101*k_[0][i]
                                                     + a103*k_[2][i]
                                                     + a104*k_[3][i]
                                                     + a105*k_[4][i]
                                                     + a106*k_[5][i]
                                                     + a107*k_[6][i]);
    ode_fun_(soltemp_, k_[9]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) {
        solemb_[i] = sol_[i] + dt*(d1*k_[0][i]
                                 + d4*k_[3][i]
                                 + d5*k_[4][i]
                                 + d6*k_[5][i]
                                 + d7*k_[6][i]
                                 + d10*k_[9][i]);
        sol_[i] = sol_[i] + dt*(b1*k_[0][i]
                              + b4*k_[3][i]
                              + b5*k_[4][i]
                              + b6*k_[5][i]
                              + b7*k_[6][i]
                              + b8*k_[7][i]
                              + b9*k_[8][i]);
    }
}
