import numpy as np

def _isclose(a, b, tol=1e-11):
    #check if numbers are close, as measured by their relative difference
    if a == b:
        return(True)
    if abs((a - b)/((a + b)/2)) < tol:
        return(True)
    return(False)

def _min2(a, b):
    #return the lesser of two numbers
    if a < b: return(a)
    return(b)

class OdeError(Exception):
    #error class for ode specific issues
    def __init__(self, message): self.message = message
    def __str__(self): return(self.message)

class OdeBase:
    """This is a base class storing information common to all the ODE
    solvers, like the number of equations in the system, the independent
    variable, counters for steps and function evaluations, etc."""

    def __init__(self, odefunk, solinit):
        #number of equations in system
        self.neq = len(solinit)
        #time (or other independent variable)
        self.t = 0.0
        #time step
        self.dt = np.nan
        #number of steps
        self.nstep = 0
        #number of function evaluations
        self.neval = 0
        #solution vector
        self.sol = np.zeros((self.neq,))
        #system of equations
        self.odefunk = odefunk
        #initial conditions
        self.solinit = np.array(solinit)
        self.sol[:] = self.solinit

class OdeBaseRK(OdeBase):
    """This class implements the non-adaptive solver functions common to
    single-step methods with three different options for output."""

    def __init__(self, odefunk, solinit):
        """
        args:
            odefunk - a function receiving the independent variable, the
                      current solution vector, and a vector to write output
                      into. This function defines the system of ODEs. Example:

                        odefunk(time, y, dydt, *args)

            solinit - a vector of initial values. The length of this vector is
                      taken as the number of equations in the system and odefunk
                      should operate on the same number of equations"""

        OdeBase.__init__(self, odefunk, solinit)
        #storage array
        self._soltemp = np.empty((self.neq,))

    def _snap(self, tout, solout, isnap, tsnap, dtsnap):
        #store a snapshot and update counters and trackers
        tout[isnap] = self.t
        solout[isnap,:] = self.sol
        return(isnap + 1, tsnap + dtsnap)

    def _solve_done(self, dt, tend):
        #test whether a solution is finished, trickier than you might think
        if self.t + dt < tend and not _isclose(self.t + dt, tend): return(False)
        return(True)

    def solve_fixed(self, tend, dt, args=(), snaps=None):
        """Integrate the system with a fixed time step
        args:
            tend - the final value of the independent variable
            dt - the size of the step to use
        optional args:
            args - additional arguments passed to odefunk
            snaps - determines the output mode
                        if snaps == 'all', all time steps are stored and
                            returned along with independent variable values
                        if snaps ==  int , evenly spaced number of snapshots
                            are stored and returned
                        if snaps == None, no output is returned"""

        #check timing
        assert tend > self.t, "tend must be greater than t, can't integrate backward"
        #check extra arguments
        assert type(args) is tuple, "args must be a tuple of extra arguments for odefunk"

        #store the time step
        self.dt = dt

        #run the solver with one of the output options
        if snaps is None:

            while self.t + self.dt < tend and not _isclose(self.t + self.dt, tend):
                self.step(self.dt, args)
            self.step(tend - self.t, args)

            return(float(self.t), self.sol.copy())

        elif snaps == 'all':

            n = int(np.ceil((tend - self.t)/self.dt)) + 1
            solout = np.zeros((n, self.neq))
            tout = np.zeros((n,))
            i = 0

            tout[i] = self.t
            solout[i,:] = self.sol
            while self.t + self.dt < tend and not _isclose(self.t + self.dt, tend):
                self.step(self.dt, args)
                i += 1
                tout[i] = self.t
                solout[i,:] = self.sol
            self.step(tend - self.t, args)
            i += 1
            tout[i] = self.t
            solout[i,:] = self.sol

            return(tout, solout)

        elif type(snaps) is int:

            tsnap = 0.0
            dtsnap = tend/float(snaps - 1)
            isnap = 0

            solout = np.zeros((snaps, self.neq))
            tout = np.zeros((snaps,))
            isnap, tsnap = self._snap(tout, solout, isnap, tsnap, dtsnap)
            while not self._solve_done(self.dt, tend) or not _isclose(tsnap, tend):

                if tsnap - self.t < self.dt or _isclose(tsnap, self.t + self.dt):
                    self.step(tsnap - self.t, args)
                    isnap, tsnap = self._snap(tout, solout, isnap, tsnap, dtsnap)
                else:
                    self.step(self.dt, args)

            self.step(tend - self.t, args)
            self._snap(tout, solout, isnap, tsnap, dtsnap)

            return(tout, solout)

        else:
            raise OdeError('the snaps argument must be a positive integer, "all", or None')

class OdeBaseARK(OdeBaseRK):
    """This class builds on the fixed step solvers in OdeBaseRK and implements
    methods for adaptive time stepping"""

    def __init__(self, odefunk, solinit, errord):
        """
        args:
            odefunk - a function receiving the independent variable, the
                      current solution vector, and a vector to write output
                      into. This function defines the system of ODEs. Example:

                        odefunk(time, y, dydt, *args)

            solinit - a vector of initial values. The length of this vector is
                      taken as the number of equations in the system and odefunk
                      should operate on the same number of equations
            errord - order of the LOWER order embedded solution"""
        OdeBaseRK.__init__(self, odefunk, solinit)
        #order of the LOWER order solution
        self.errord = errord
        #low order solution
        self._sollo = np.empty((self.neq,))
        #high order solution
        self._solhi = np.empty((self.neq,))
        #storage for cancelling a step and reversing
        self._solprev = np.empty((self.neq,))
        #maximum allowable time step
        self.dtmax = np.inf
        #factors for adaptive time stepping
        self.facsafe = 0.9
        self.facmin = 1.0/100
        self.facmax = 10.0
        self.abstol = 1e-6
        self.reltol = 1e-6
        self.nrej = 0

    def _error_est(self):
        d = np.abs(self._sollo - self._solhi)
        sc = self.abstol + self.reltol*np.maximum(np.abs(self._sollo), np.abs(self._solhi))
        err = np.sqrt(((d*d)/(sc*sc)).sum())
        return(err)

    def _facopt(self, err):
        if err == 0: return(self.facmax)
        facopt = self.facsafe*(1./err)**(1./(self.errord + 1.0))
        if facopt > self.facmax: return(self.facmax)
        if facopt < self.facmin: return(self.facmin)
        return(facopt)

    def _reverse_step(self):
        np.copyto(self.sol, self._solprev)
        self.t -= self.dt
        self.nrej += 1
        self.nstep -= 1

    def solve_adaptive(self, tend, dt0, args=(), snaps=None):
        """Integrate the system with an adaptive time step. The adapting
        variables (abstol, reltol, facmax, ...) are class variables.
        args:
            tend - the final value of the independent variable
            dt0 - the initial time step to use
        optional args:
            args - additional arguments passed to odefunk
            snaps - determines the output mode
                        if snaps == 'all', all time steps are stored and
                            returned along with independent variable values
                        if snaps ==  int , evenly spaced number of snapshots
                            are stored and returned
                        if snaps == None, no output is returned"""

        #check timing
        assert tend >= self.t, "tend must be greater than or equal to t, can't integrate backward"
        #check extra arguments
        assert type(args) is tuple, "args must be a tuple of extra arguments for odefunk"

        #store the time step
        self.dt = dt0

        #run the solver with one of the output options
        if snaps is None:

            while self.t < tend:
                self.step_adaptive(self.dt, args)
                err = self._error_est()
                if err > 1.0:
                    self._reverse_step()
                #compute optimal time step, can't go past tend
                self.dt = _min2(self.dt*self._facopt(err), tend - self.t)
                #can't be larger than dtmax
                self.dt = _min2(self.dt, self.dtmax)
            return(float(self.t), self.sol.copy())

        elif snaps == 'all':

            tout = [self.t]
            solout = [self.sol.copy()]
            while self.t < tend:
                self.step_adaptive(self.dt, args)
                err = self._error_est()
                if err > 1.0:
                    self._reverse_step()
                else:
                    tout.append(self.t)
                    solout.append(self.sol.copy())
                #compute optimal time step, can't go past tend
                self.dt = _min2(self.dt*self._facopt(err), tend - self.t)
                #can't be larger than dtmax
                self.dt = _min2(self.dt, self.dtmax)
            return(tout, np.vstack(solout))

        elif type(snaps) is int:

            #initialize snap trackers
            tsnap = 0.0
            dtsnap = tend/float(snaps - 1)
            isnap = 0
            #initialize solution variables
            solout = np.zeros((snaps, self.neq))
            tout = np.zeros((snaps,))
            isnap, tsnap = self._snap(tout, solout, isnap, tsnap, dtsnap)
            #solve
            while self.t < tend:
                self.step_adaptive(self.dt, args)
                err = self._error_est()
                if err > 1.0:
                    self._reverse_step()
                elif _isclose(self.t, tsnap):
                    isnap, tsnap = self._snap(tout, solout, isnap, tsnap, dtsnap)
                #compute optimal time step, can't go past tend
                self.dt = _min2(self.dt*self._facopt(err), tend - self.t)
                #can't be larger than dtmax
                self.dt = _min2(self.dt, self.dtmax)
                #can't go past the next snap point
                self.dt = _min2(self.dt, tsnap - self.t)
            return(tout, solout)

        else:
            raise OdeError('the snaps argument must be a positive integer or "all"')

class OdeBaseIRK(OdeBaseRK):

    def __init__(self, odefunk, solinit):
        OdeBase.__init__(self, odefunk, solinit)
