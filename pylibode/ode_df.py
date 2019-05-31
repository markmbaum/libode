import pandas as pd

#import all the solver classes
from .ode_explicit_RK import *
from .ode_implicit_RK import *

#put the solvers together with some other information
solvers = [('Euler', "Euler's Method", OdeEuler, 1, False, 'explicit'),
           ('Trapz', 'Trapezoidal Method', OdeTrapz, 2, True, 'explicit'),
           ('Midpnt', 'Midpoint Method', OdeMidpnt, 2, True, 'explicit'),
           ('Ssp3', 'Strong Stability Preserving Method of Order 3', OdeSsp3, 3, False, 'explicit'),
           ('RKF32', 'Runge-Kutta Fehlberg 3(2) Method', OdeRKF32, 3, True, 'explicit'),
           ('RK4', 'Classic Runge-Kutta 4th Order Method', OdeRK4, 4, False, 'explicit'),
           ('DoPri54', 'Dormand-Prince 5(4) Method', OdeDoPri54, 5, True, 'explicit'),
           ('Butcher6', "Butcher's 6th Order Method", OdeButcher6, 6, False, 'explicit'),
           ('DoPri87', 'Dormand-Prince 8(7) Method', OdeDoPri87, 8, True, 'explicit'),
           ('BackwardEuler', "Implicit Euler's Method (Backward Euler)", OdeBackwardEuler, 1, False, 'implicit'),
           ('ImplicitMidpnt', 'Implicit Midpoint Method', OdeImplicitMidpnt, 2, False, 'implicit'),
           ('ImplicitTrapz', 'Implicit Trapezoidal Method', OdeImplicitTrapz, 2, False, 'implicit')]

#package them up into a dataframe
solvers = pd.DataFrame(
    dict(
        zip(
            ['name', 'long_name', 'constructor', 'order', 'adaptive', 'type'],
            zip(*solvers)
        )
    )
)
solvers.set_index('name', inplace=True)
