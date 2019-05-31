import numpy as np

from .ode_bases import OdeBaseRK, OdeBaseARK

class OdeEuler(OdeBaseRK):
    """Euler's Method, 1st order, explicit, non-adaptive
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseRK.__init__(self, odefunk, solinit)
        #k vectors
        self._k = np.empty((self.neq,))

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k, *args)
        self.neval += 1

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*self._k

        self.t += dt
        self.nstep += 1

class OdeTrapz(OdeBaseARK):
    """Trapezoidal Method, 2nd order, explicit, adaptive
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseARK.__init__(self, odefunk, solinit, 1)
        #k vectors
        self._k1 = np.empty((self.neq,))
        self._k2 = np.empty((self.neq,))
        #tableau of coefficients
        self._c2 = 1.0; self._a21 = 1.0
        self._b1 = 1.0/2; self._b2 = 1.0/2

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k1, *args)
        self._sol_temp = self.sol + dt*self._a21*self._k1

        self.odefunk(self.t + dt*self._c2, self._sol_temp, self._k2, *args)

        self.neval += 2

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*(self._b1*self._k1 + self._b2*self._k2)

        self.t += dt
        self.nstep += 1

    def step_adaptive(self, dt, args):

        self._calc_k(dt, args)

        self._solhi = self.sol + dt*(self._b1*self._k1 + self._b2*self._k2)
        self._sollo = self.sol + dt*self._k1

        np.copyto(self._solprev, self.sol)
        np.copyto(self.sol, self._solhi)

        self.t += dt
        self.nstep += 1

class OdeMidpnt(OdeBaseARK):
    """Midpoint Method, 2nd order, explicit, adaptive
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseARK.__init__(self, odefunk, solinit, 1)
        #k vectors
        self._k1 = np.empty((self.neq,))
        self._k2 = np.empty((self.neq,))
        #tableau of coefficients
        self._c2 = 1.0/2; self._a21 = 1.0/2
        self._b1 = 0.0; self._b2 = 1.0

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k1, *args)
        self._sol_temp = self.sol + dt*self._a21*self._k1

        self.odefunk(self.t + dt*self._c2, self._sol_temp, self._k2, *args)

        self.neval += 2

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*self._b2*self._k2

        self.t += dt
        self.nstep += 1

    def step_adaptive(self, dt, args):

        self._calc_k(dt, args)

        self._solhi = self.sol + dt*self._b2*self._k2
        self._sollo = self.sol + dt*self._k1

        np.copyto(self._solprev, self.sol)
        np.copyto(self.sol, self._solhi)

        self.t += dt
        self.nstep += 1

class OdeSsp3(OdeBaseRK):
    """Strong Stability Preserving Method, 3rd order, explicit, non-adaptive
    C. W. Shu and S. Osher, Effcient implementation of essentially nonoscillatory shock-capturing schemes, J. Comput. Phys., 77, 1988, pp. 439-471.
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseRK.__init__(self, odefunk, solinit)
        #k vectors
        self._k1 = np.empty((self.neq,))
        self._k2 = np.empty((self.neq,))
        self._k3 = np.empty((self.neq,))
        #tableau of coefficients
        self._c2 = 1.0; self._a21 = 1.0
        self._c3 = 1.0/2; self._a31 = 1.0/4; self._a32 = 1.0/4
        self._b1 = 1.0/6; self._b2 = 1.0/6; self._b3 = 2.0/3

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k1, *args)
        self._sol_temp = self.sol + dt*self._a21*self._k1

        self.odefunk(self.t + dt*self._c2, self._sol_temp, self._k2, *args)
        self._sol_temp = self.sol + dt*(self._a31*self._k1 + self._a32*self._k2)

        self.odefunk(self.t + dt*self._c3, self._sol_temp, self._k3, *args)

        self.neval += 3

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*(self._b1*self._k1 + self._b2*self._k2 + self._b3*self._k3)

        self.t += dt
        self.nstep += 1

class OdeRKF32(OdeBaseARK):
    """Runge-Kutta-Fehlberg Method, 3rd order, explicit, adaptive
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseARK.__init__(self, odefunk, solinit, 2)
        #k vectors
        self._k1 = np.empty((self.neq,))
        self._k2 = np.empty((self.neq,))
        self._k3 = np.empty((self.neq,))
        #tableau of coefficients
        self._c2 = 1.0; self._a21 = 1.0
        self._c3 = 1.0/2; self._a31 = 1.0/4; self._a32 = 1.0/4
        self._b1 = 1.0/2; self._b2  = 1.0/2; self._b3 = 0.0
        self._d1 = 1.0/6; self._d2 = 1.0/6; self._d3  = 4.0/6

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k1, *args)

        self._soltemp = self.sol + dt*self._a21*self._k1
        self.odefunk(self.t + dt*self._c2, self._soltemp, self._k2, *args)

        self._soltemp = self.sol + dt*(self._a31*self._k1 + self._a32*self._k2)
        self.odefunk(self.t + dt*self._c3, self._soltemp, self._k3, *args)

        self.neval += 3

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*(self._d1*self._k1 + self._d2*self._k2 + self._d3*self._k3)

        self.t += dt
        self.nstep += 1

    def step_adaptive(self, dt, args):

        self._calc_k(dt, args)

        self._sollo = self.sol + dt*(self._b1*self._k1 + self._b2*self._k2 + self._b3*self._k3)
        self._solhi = self.sol + dt*(self._d1*self._k1 + self._d2*self._k2 + self._d3*self._k3)

        np.copyto(self._solprev, self.sol)
        np.copyto(self.sol, self._solhi)

        self.t += dt
        self.nstep += 1

class OdeRK4(OdeBaseRK):
    """Classic RK, 4th order, explicit, non-adaptive
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseRK.__init__(self, odefunk, solinit)
        #k vectors
        self._k1 = np.empty((self.neq,))
        self._k2 = np.empty((self.neq,))
        self._k3 = np.empty((self.neq,))
        self._k4 = np.empty((self.neq,))
        #tableau of coefficients
        self._c2 = 1.0/2; self._a21 = 1.0/2
        self._c3 = 1.0/2; self._a32 = 1.0/2
        self._c4 = 1.0; self._a43 = 1.0
        self._b1 = 1.0/6; self._b2 = 1.0/3; self._b3 = 1.0/3; self._b4 = 1.0/6

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k1, *args)
        self._sol_temp = self.sol + dt*self._a21*self._k1

        self.odefunk(self.t + dt*self._c2, self._sol_temp, self._k2, *args)
        self._sol_temp = self.sol + dt*self._a32*self._k2

        self.odefunk(self.t + dt*self._c3, self._sol_temp, self._k3, *args)
        self._sol_temp = self.sol + dt*self._a43*self._k3

        self.odefunk(self.t + dt*self._c4, self._sol_temp, self._k4, *args)

        self.neval += 4

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*(self._b1*self._k1 + self._b2*self._k2
                                + self._b3*self._k3 + self._b4*self._k4)

        self.t += dt
        self.nstep += 1

class OdeDoPri54(OdeBaseARK):
    """Dormand-Prince Method, 5th order, explicit, adaptive
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseARK.__init__(self, odefunk, solinit, 4)
        #k vectors
        self._k1 = np.empty((self.neq,))
        self._k2 = np.empty((self.neq,))
        self._k3 = np.empty((self.neq,))
        self._k4 = np.empty((self.neq,))
        self._k5 = np.empty((self.neq,))
        self._k6 = np.empty((self.neq,))
        self._k7 = np.empty((self.neq,))
        #tableau of coefficents
        self._c2 = 1.0/5
        self._a21 = 1.0/5
        self._c3 = 3.0/10
        self._a31 = 3.0/40
        self._a32 = 9.0/40
        self._c4 = 4.0/5
        self._a41 = 44.0/45
        self._a42 = -56.0/15
        self._a43 = 32.0/9
        self._c5 = 8.0/9
        self._a51 = 19372.0/6561
        self._a52 = -25360.0/2187
        self._a53 = 64448.0/6561
        self._a54 = -212.0/729
        self._c6 = 1.0
        self._a61 = 9017.0/3168
        self._a62 = -355.0/33
        self._a63 = 46732.0/5247
        self._a64 = 49.0/176
        self._a65 = -5103.0/18656
        self._c7 = 1.0
        self._a71 = 35.0/384
        self._a72 = 0.0
        self._a73 = 500.0/1113
        self._a74 = 125.0/192
        self._a75 = -2187.0/6784
        self._a76 = 11.0/84
        self._b1 = 35.0/384
        self._b2  = 0.0
        self._b3  = 500.0/1113
        self._b4  = 125.0/192
        self._b5  = -2187.0/6784
        self._b6 = 11.0/84
        self._d1 = 5179.0/57600
        self._d2  = 0.0
        self._d3  = 7571.0/16695
        self._d4 = 393.0/640
        self._d5  = -92097.0/339200
        self._d6 = 187.0/2100
        self._d7 = 1.0/40

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k1, *args)

        self._soltemp = self.sol + dt*self._a21*self._k1
        self.odefunk(self.t + dt*self._c2, self._soltemp, self._k2, *args)

        self._soltemp = self.sol + dt*(self._a31*self._k1 + self._a32*self._k2)
        self.odefunk(self.t + dt*self._c3, self._soltemp, self._k3, *args)

        self._soltemp = self.sol + dt*(self._a41*self._k1 + self._a42*self._k2 + self._a43*self._k3)
        self.odefunk(self.t + dt*self._c4, self._soltemp, self._k4, *args)

        self._soltemp = self.sol + dt*(self._a51*self._k1 + self._a52*self._k2
                                    + self._a53*self._k3 + self._a54*self._k4)
        self.odefunk(self.t + dt*self._c5, self._soltemp, self._k5, *args)

        self._soltemp = self.sol + dt*(self._a61*self._k1 + self._a62*self._k2
                                    + self._a63*self._k3 + self._a64*self._k4
                                    + self._a65*self._k5)
        self.odefunk(self.t + dt*self._c6, self._soltemp, self._k6, *args)

        self._soltemp = self.sol + dt*(self._a71*self._k1 + self._a72*self._k2
                                    + self._a73*self._k3 + self._a74*self._k4
                                    + self._a75*self._k5 + self._a76*self._k6)
        self.odefunk(self.t + dt*self._c7, self._soltemp, self._k7, *args)

        self.neval += 6

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*(self._b1*self._k1 + self._b3*self._k3
                                    + self._b4*self._k4 + self._b5*self._k5
                                    + self._b6*self._k6)

        self.t += dt
        self.nstep += 1

    def step_adaptive(self, dt, args):

        self._calc_k(dt, args)

        self._solhi = self.sol + dt*(self._b1*self._k1 + self._b3*self._k3
                                    + self._b4*self._k4 + self._b5*self._k5
                                    + self._b6*self._k6)
        self._sollo = self.sol + dt*(self._d1*self._k1 + self._d3*self._k3
                                    + self._d4*self._k4 + self._d5*self._k5
                                    + self._d6*self._k6 + self._d7*self._k7)

        np.copyto(self._solprev, self.sol)
        np.copyto(self.sol, self._solhi)

        self.t += dt
        self.nstep += 1

class OdeButcher6(OdeBaseRK):
    """Butcher's Method, 6th order, explicit, non-adaptive
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseRK.__init__(self, odefunk, solinit)
        #k vectors
        self._k1 = np.empty((self.neq,))
        self._k2 = np.empty((self.neq,))
        self._k3 = np.empty((self.neq,))
        self._k4 = np.empty((self.neq,))
        self._k5 = np.empty((self.neq,))
        self._k6 = np.empty((self.neq,))
        self._k7 = np.empty((self.neq,))
        #tableau of coefficents
        self._c2 = 1.0/2
        self._a21 = 1.0/2
        self._c3 = 2.0/3
        self._a31 = 2.0/9
        self._a32 = 4.0/9
        self._c4 = 1.0/3
        self._a41 = 7.0/36
        self._a42 = 2.0/9
        self._a43 = -1.0/12
        self._c5 = 5.0/6
        self._a51 = -35.0/144
        self._a52 = -55.0/36
        self._a53 = 35.0/48
        self._a54 = 15.0/8
        self._c6 = 1.0/6
        self._a61 = -1.0/360
        self._a62 = -11.0/36
        self._a63 = -1.0/8
        self._a64 = 1.0/2
        self._a65 = 1.0/10
        self._c7 = 1.0
        self._a71 = -41.0/260
        self._a72 = 22.0/13
        self._a73 = 43.0/156
        self._a74 = -118.0/39
        self._a75 = 32.0/195
        self._a76 = 80.0/39
        self._b1 = 13.0/200
        self._b2 = 0.0
        self._b3 = 11.0/40
        self._b4 = 11.0/40
        self._b5 = 4.0/25
        self._b6 = 4.0/25
        self._b7 = 13.0/200

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k1, *args)

        self._soltemp = self.sol + dt*self._a21*self._k1
        self.odefunk(self.t + dt*self._c2, self._soltemp, self._k2, *args)

        self._soltemp = self.sol + dt*(self._a31*self._k1 + self._a32*self._k2)
        self.odefunk(self.t + dt*self._c3, self._soltemp, self._k3, *args)

        self._soltemp = self.sol + dt*(self._a41*self._k1 + self._a42*self._k2
                                    + self._a43*self._k3)
        self.odefunk(self.t + dt*self._c4, self._soltemp, self._k4, *args)

        self._soltemp = self.sol + dt*(self._a51*self._k1 + self._a52*self._k2
                                    + self._a53*self._k3 + self._a54*self._k4)
        self.odefunk(self.t + dt*self._c5, self._soltemp, self._k5, *args)

        self._soltemp = self.sol + dt*(self._a61*self._k1 + self._a62*self._k2
                                    + self._a63*self._k3 + self._a64*self._k4
                                    + self._a65*self._k5)
        self.odefunk(self.t + dt*self._c6, self._soltemp, self._k6, *args)

        self._soltemp = self.sol + dt*(self._a71*self._k1 + self._a72*self._k2
                                    + self._a73*self._k3 + self._a74*self._k4
                                    + self._a75*self._k5 + self._a76*self._k6)
        self.odefunk(self.t + dt*self._c7, self._soltemp, self._k7, *args)

        self.neval += 7

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*(self._b1*self._k1 + self._b3*self._k3
                                    + self._b4*self._k4 + self._b5*self._k5
                                    + self._b6*self._k6 + self._b7*self._k7)

        self.t += dt
        self.nstep += 1

class OdeDoPri87(OdeBaseARK):
    """Dormand-Prince Method, 8th order, explicit, adaptive
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseARK.__init__(self, odefunk, solinit, 7)
        #k vectors
        self._k1 = np.empty((self.neq,))
        self._k2 = np.empty((self.neq,))
        self._k3 = np.empty((self.neq,))
        self._k4 = np.empty((self.neq,))
        self._k5 = np.empty((self.neq,))
        self._k6 = np.empty((self.neq,))
        self._k7 = np.empty((self.neq,))
        self._k8 = np.empty((self.neq,))
        self._k9 = np.empty((self.neq,))
        self._k10 = np.empty((self.neq,))
        self._k11 = np.empty((self.neq,))
        self._k12 = np.empty((self.neq,))
        self._k13 = np.empty((self.neq,))
        #tableau of coefficients
        self._c2 = 1.0/18
        self._a21 = 1.0/18
        self._c3 = 1.0/12
        self._a31 = 1.0/48
        self._a32 = 1.0/16
        self._c4 = 1.0/8
        self._a41 = 1.0/32
        self._a42 = 0.0
        self._a43 = 3.0/32
        self._c5 = 5.0/16
        self._a51 = 5.0/16
        self._a52 = 0.0
        self._a53 = -75.0/64
        self._a54 = 75.0/64
        self._c6 = 3.0/8
        self._a61 = 3.0/80
        self._a62 = 0.0
        self._a63 = 0.0
        self._a64 = 3.0/16
        self._a65 = 3.0/20
        self._c7 = 59.0/400
        self._a71 = 29443841.0/614563906
        self._a72 = 0.0
        self._a73 = 0.0
        self._a74 = 77736538.0/692538347
        self._a75 = -28693883.0/1125000000
        self._a76 = 23124283.0/1800000000
        self._c8 = 93.0/200
        self._a81 = 16016141.0/946692911
        self._a82 = 0.0
        self._a83 = 0.0
        self._a84 = 61564180.0/158732637
        self._a85 = 22789713.0/633445777
        self._a86 = 545815736.0/2771057229
        self._a87 = -180193667.0/1043307555
        self._c9 = 5490023248.0/9719169821
        self._a91 = 39632708.0/573591083
        self._a92 = 0.0
        self._a93 = 0.0
        self._a94 = -433636366.0/683701615
        self._a95 = -421739975.0/2616292301
        self._a96 = 100302831.0/723423059
        self._a97 = 790204164.0/839813087
        self._a98 = 800635310.0/3783071287
        self._c10 = 13.0/20
        self._a101 = 246121993.0/1340847787
        self._a102 = 0.0
        self._a103 = 0.0
        self._a104 = -37695042795.0/15268766246
        self._a105 = -309121744.0/1061227803
        self._a106 = -12992083.0/490766935
        self._a107 = 6005943493.0/2108947869
        self._a108 = 393006217.0/1396673457
        self._a109 = 123872331.0/1001029789
        self._c11 = 1201146811.0/1299019798
        self._a111 = -1028468189.0/846180014
        self._a112 = 0.0
        self._a113 = 0.0
        self._a114 = 8478235783.0/508512852
        self._a115 = 1311729495.0/1432422823
        self._a116 = -10304129995.0/1701304382
        self._a117 = -48777925059.0/3047939560
        self._a118 = 15336726248.0/1032824649
        self._a119 = -45442868181.0/3398467696
        self._a1110 = 3065993473.0/597172653
        self._c12 = 1.0
        self._a121 = 185892177.0/718116043
        self._a122 = 0.0
        self._a123 = 0.0
        self._a124 = -3185094517.0/667107341
        self._a125 = -477755414.0/1098053517
        self._a126 = -703635378.0/230739211
        self._a127 = 5731566787.0/1027545527
        self._a128 = 5232866602.0/850066563
        self._a129 = -4093664535.0/808688257
        self._a1210 = 3962137247.0/1805957418
        self._a1211 = 65686358.0/487910083
        self._c13 = 1.0
        self._a131 = 403863854.0/491063109
        self._a132 = 0.0
        self._a133 = 0.0
        self._a134 = -5068492393.0/434740067
        self._a135 = -411421997.0/543043805
        self._a136 = 652783627.0/914296604
        self._a137 = 11173962825.0/925320556
        self._a138 = -13158990841.0/6184727034
        self._a139 = 3936647629.0/1978049680
        self._a1310 = -160528059.0/685178525
        self._a1311 = 248638103.0/1413531060
        self._a1312 = 0.0
        self._b1 = 14005451.0/335480064
        self._b2 = 0.0
        self._b3 = 0.0
        self._b4 = 0.0
        self._b5 = 0.0
        self._b6 = -59238493.0/1068277825
        self._b7 = 181606767.0/758867731
        self._b8 = 561292985.0/797845732
        self._b9 = -1041891430.0/1371343529
        self._b10 = 760417239.0/1151165299
        self._b11 = 118820643.0/751138087
        self._b12 = -528747749.0/2220607170
        self._b13 = 1.0/4
        self._d1 = 13451932.0/455176632
        self._d2 = 0.0
        self._d3 = 0.0
        self._d4 = 0.0
        self._d5 = 0.0
        self._d6 = -808719846.0/976000145
        self._d7 = 1757004468.0/5645159321
        self._d8 = 656045339.0/265891186
        self._d9 = -3867574721.0/1518517206
        self._d10 = 465885868.0/322736535
        self._d11 = 53011238.0/667516719
        self._d12 = 2.0/45
        self._d13 = 0.0

    def _calc_k(self, dt, args):

        self.odefunk(self.t, self.sol, self._k1, *args)

        self._soltemp = self.sol + dt*self._a21*self._k1
        self.odefunk(self.t + dt*self._c2, self._soltemp, self._k2, *args)

        self._soltemp = self.sol + dt*(self._a31*self._k1 + self._a32*self._k2)
        self.odefunk(self.t + dt*self._c3, self._soltemp, self._k3, *args)

        self._soltemp = self.sol + dt*(self._a41*self._k1 + self._a43*self._k3)
        self.odefunk(self.t + dt*self._c4, self._soltemp, self._k4, *args)

        self._soltemp = self.sol + dt*(self._a51*self._k1 + self._a53*self._k3 + self._a54*self._k4)
        self.odefunk(self.t + dt*self._c5, self._soltemp, self._k5, *args)

        self._soltemp = self.sol + dt*(self._a61*self._k1 + self._a64*self._k4 + self._a65*self._k5)
        self.odefunk(self.t + dt*self._c6, self._soltemp, self._k6, *args)

        self._soltemp = self.sol + dt*(self._a71*self._k1 + self._a74*self._k4
                                    + self._a75*self._k5 + self._a76*self._k6)
        self.odefunk(self.t + dt*self._c7, self._soltemp, self._k7, *args)

        self._soltemp = self.sol + dt*(self._a81*self._k1 + self._a84*self._k4
                                    + self._a85*self._k5 + self._a86*self._k6
                                    + self._a87*self._k7)
        self.odefunk(self.t + dt*self._c8, self._soltemp, self._k8, *args)

        self._soltemp = self.sol + dt*(self._a91*self._k1 + self._a94*self._k4
                                    + self._a95*self._k5 + self._a96*self._k6
                                    + self._a97*self._k7 + self._a98*self._k8)
        self.odefunk(self.t + dt*self._c9, self._soltemp, self._k9, *args)

        self._soltemp = self.sol + dt*(self._a101*self._k1 + self._a104*self._k4
                                    + self._a105*self._k5 + self._a106*self._k6
                                    + self._a107*self._k7 + self._a108*self._k8
                                    + self._a109*self._k9)
        self.odefunk(self.t + dt*self._c10, self._soltemp, self._k10, *args)

        self._soltemp = self.sol + dt*(self._a111*self._k1 + self._a114*self._k4
                                    + self._a115*self._k5 + self._a116*self._k6
                                    + self._a117*self._k7 + self._a118*self._k8
                                    + self._a119*self._k9 + self._a1110*self._k10)
        self.odefunk(self.t + dt*self._c11, self._soltemp, self._k11, *args)

        self._soltemp = self.sol + dt*( self._a121*self._k1 + self._a124*self._k4
                                    + self._a125*self._k5 + self._a126*self._k6
                                    + self._a127*self._k7
                                    + self._a128*self._k8 + self._a129*self._k9
                                    + self._a1210*self._k10 + self._a1211*self._k11)
        self.odefunk(self.t + dt*self._c12, self._soltemp, self._k12, *args)

        self._soltemp = self.sol + dt*( self._a131*self._k1 + self._a134*self._k4
                                    + self._a135*self._k5 + self._a136*self._k6
                                    + self._a137*self._k7
                                    + self._a138*self._k8 + self._a139*self._k9
                                    + self._a1310*self._k10 + self._a1311*self._k11)
        self.odefunk(self.t + dt*self._c13, self._soltemp, self._k13, *args)

        self.neval += 13

    def step(self, dt, args):

        self._calc_k(dt, args)

        self.sol = self.sol + dt*(self._b1*self._k1 + self._b6*self._k6
                                + self._b7*self._k7 + self._b8*self._k8
                                + self._b9*self._k9 + self._b10*self._k10
                                + self._b11*self._k11 + self._b12*self._k12
                                + self._b13*self._k13)

        self.t += dt
        self.nstep += 1

    def step_adaptive(self, dt, args):

        self._calc_k(dt, args)

        self._sollo = self.sol + dt*(self._d1*self._k1 + self._d6*self._k6
                                    + self._d7*self._k7 + self._d8*self._k8
                                    + self._d9*self._k9 + self._d10*self._k10
                                    + self._d11*self._k11 + self._d12*self._k12
                                    + self._d13*self._k13);
        self._solhi = self.sol + dt*(self._b1*self._k1 + self._b6*self._k6
                                    + self._b7*self._k7 + self._b8*self._k8
                                    + self._b9*self._k9 + self._b10*self._k10
                                    + self._b11*self._k11 + self._b12*self._k12
                                    + self._b13*self._k13)

        np.copyto(self._solprev, self.sol)
        np.copyto(self.sol, self._solhi)

        self.t += dt
        self.nstep += 1
