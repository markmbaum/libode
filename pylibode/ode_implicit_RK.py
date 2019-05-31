import numpy as np
from scipy.optimize import root

from .ode_bases import OdeBaseIRK

def _funk_BackwardEuler(x, s, dt):
    s.odefunk(s.t + dt, x, s._k)
    s.neval += 1
    s._soltemp = s._k*dt + s.sol
    return(s._soltemp - x)

class OdeBackwardEuler(OdeBaseIRK):
    """Backward Euler, 1st order, implicit
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseIRK.__init__(self, odefunk, solinit)
        self._k = np.zeros((self.neq,))

    def step(self, dt):
        self.sol[:] = root(_funk_BackwardEuler, self.sol, (self, dt)).x

        self.t += dt
        self.nstep += 1

def _funk_ImplicitMidpnt(x, s, dt):
    s.odefunk(s.t + dt/2, (s.sol + x)/2, s._k)
    s.neval += 1
    s._soltemp = s._k*dt + s.sol
    return(s._soltemp - x)

class OdeImplicitMidpnt(OdeBaseIRK):
    """Midpoint Method, 2nd order, implicit
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseIRK.__init__(self, odefunk, solinit)
        self._k = np.zeros((self.neq,))

    def step(self, dt):
        self.sol[:] = root(_funk_ImplicitMidpnt, self.sol, (self, dt)).x

        self.t += dt
        self.nstep += 1

def _funk_ImplicitTrapz(x, s, dt):
    s.odefunk(s.t + dt, x, s._k2)
    s.neval += 1
    s._soltemp = (s._k1 + s._k2)*dt/2 + s.sol
    return(s._soltemp - x)

class OdeImplicitTrapz(OdeBaseIRK):
    """Implicit Trapezoidal Method, 2nd order, implicit
    init args:
        odefunk - a function receiving the independent variable, the current
                  solution vector, and a vector to write output into. This
                  function defines the system of ODEs. For example:

                    odefunk(time, y, dydt, *args)

        solinit - a vector of initial values. The length of this vector is
                  taken as the number of equations in the system and odefunk
                  should operate on the same number of equations"""

    def __init__(self, odefunk, solinit):
        OdeBaseIRK.__init__(self, odefunk, solinit)
        self._k1 = np.zeros((self.neq,))
        self._k2 = np.zeros((self.neq,))

    def step(self, dt):
        self.odefunk(self.t, self.sol, self._k1)
        self.sol[:] = root(_funk_ImplicitTrapz, self.sol, (self, dt)).x

        self.t += dt
        self.nstep += 1
