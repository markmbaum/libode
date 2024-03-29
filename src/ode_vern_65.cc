//! \file ode_vern_65.cc

#include "ode_vern_65.h"

namespace ode {

OdeVern65::OdeVern65 (unsigned long neq) :
    OdeEmbedded (neq, false, 5),
    OdeRK (neq, 9),
    OdeERK (neq) {

    method_ = "Vern65";

    /*
    coefficients copied from:
        http://people.math.sfu.ca/~jverner/
        http://people.math.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.CoeffsOnlyFLOAT
    */

    c2 = 0.06;
    c3 = 0.09593333333333333333333333333333333333333;
    c4 = 0.1439;
    c5 = 0.4973;
    c6 = 0.9725;
    c7 = 0.9995;
    c8 = 1.0;
    c9 = 1.0;

    a21 =  c2;

    a31 =  0.01923996296296296296296296296296296296296;
    a32 =  0.07669337037037037037037037037037037037037;

    a41 =  0.035975;
    a43 =  0.107925;

    a51 =  1.318683415233148260919747276431735612861;
    a53 = -5.042058063628562225427761634715637693344;
    a54 =  4.220674648395413964508014358283902080483;

    a61 = -41.87259166432751461803757780644346812905;
    a63 =  159.4325621631374917700365669070346830453;
    a64 = -122.1192135650100309202516203389242140663;
    a65 =  5.531743066200053768252631238332999150076;

    a71 = -54.43015693531650433250642051294142461271;
    a73 =  207.0672513650184644273657173866509835987;
    a74 = -158.6108137845899991828742424365058599469;
    a75 =  6.991816585950242321992597280791793907096;
    a76 = -0.01859723106220323397765171799549294623692;

    a81 = -54.66374178728197680241215648050386959351;
    a83 =  207.9528062553893734515824816699834244238;
    a84 = -159.2889574744995071508959805871426654216;
    a85 =  7.018743740796944434698170760964252490817;
    a86 = -0.01833878590504572306472782005141738268361;
    a87 = -0.0005119484997882099077875432497245168395840;

    a91 =  0.03438957868357036009278820124728322386520;
    a94 =  0.2582624555633503404659558098586120858767;
    a95 =  0.4209371189673537150642551514069801967032;
    a96 =  4.405396469669310170148836816197095664891;
    a97 = -176.4831190242986576151740942499002125029;
    a98 =  172.3641334014150730294022582711902413315;

    b1 =  0.03438957868357036009278820124728322386520;
    b4 =  0.2582624555633503404659558098586120858767;
    b5 =  0.4209371189673537150642551514069801967032;
    b6 =  4.405396469669310170148836816197095664891;
    b7 = -176.4831190242986576151740942499002125029;
    b8 =  172.3641334014150730294022582711902413315;

    d1 =  0.04909967648382489730906854927971225836479;
    d4 =  0.2251112229516524153401395320539875329485;
    d5 =  0.4694682253029562039431948525047387412553;
    d6 =  0.8065792249988867707634161808995217981443;
    d8 = -0.6071194891777959797672951465256217122488;
    d9 =  0.05686113944047569241147603178766138153594;
}

void OdeVern65::step_ (double dt) {

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
                                                     + a94*k_[3][i]
                                                     + a95*k_[4][i]
                                                     + a96*k_[5][i]
                                                     + a97*k_[6][i]
                                                     + a98*k_[7][i]);
    ode_fun_(soltemp_, k_[8]);

    //------------------------------------------------------------------
    for (i=0; i<neq_; i++) {
        solemb_[i] = sol_[i] + dt*(d1*k_[0][i]
                                 + d4*k_[3][i]
                                 + d5*k_[4][i]
                                 + d6*k_[5][i]
                                 + d8*k_[7][i]
                                 + d9*k_[8][i]);
        sol_[i] = sol_[i] + dt*(b1*k_[0][i]
                              + b4*k_[3][i]
                              + b5*k_[4][i]
                              + b6*k_[5][i]
                              + b7*k_[6][i]
                              + b8*k_[7][i]);
    }
}

} // namespace ode
