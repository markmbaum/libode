/*
This class implements the 7th and 8th order method of Dormand and Prince, as reprinted in [1]. The tableau coefficients are terribly unweildy, but the solver is scary accurate.

    [1] E. Hairer, S. P. Nørsett, and G. Wanner. 1993. Solving Ordinary Differential Equations I (2nd Revised. Ed.): Nonstiff Problems. Springer-Verlag, Berlin, Heidelberg.
*/

#ifndef ODE_DOPRI87_HPP
#define ODE_DOPRI87_HPP

#include "OdeBaseARK.hpp"

class OdeDoPri87 : public OdeBaseARK {
public:
    //constructor
    OdeDoPri87 (unsigned long long neq_);
    //destructor
    ~OdeDoPri87 ();
    //function for taking a single time step
    void step (double dt);
private:
    //arrays for time stepping routine
    double *k1, *k2, *k3, *k4, *k5, *k6, *k7, *k8, *k9, *k10, *k11, *k12, *k13;
    //coefficents of tableau
    double c2,   a21,
           c3,   a31,   a32,
           c4,   a41,   a42,  a43,
           c5,   a51,   a52,  a53,  a54,
           c6,   a61,   a62,  a63,  a64,  a65,
           c7,   a71,   a72,  a73,  a74,  a75,  a76,
           c8,   a81,   a82,  a83,  a84,  a85,  a86,  a87,
           c9,   a91,   a92,  a93,  a94,  a95,  a96,  a97,  a98,
           c10,  a101, a102, a103, a104, a105, a106, a107, a108, a109,
           c11,  a111, a112, a113, a114, a115, a116, a117, a118, a119, a1110,
           c12,  a121, a122, a123, a124, a125, a126, a127, a128, a129, a1210, a1211,
           c13,  a131, a132, a133, a134, a135, a136, a137, a138, a139, a1310, a1311, a1312,
                   b1,   b2,   b3,   b4,   b5,   b6,   b7,   b8,   b9,   b10,   b11,   b12,   b13,
                   d1,   d2,   d3,   d4,   d5,   d6,   d7,   d8,   d9,   d10,   d11,   d12,   d13;

};

#endif
