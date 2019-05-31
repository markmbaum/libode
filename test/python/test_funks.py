import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from numpy import cos, sin, log, exp

#-------------------------------------------------------------------------------
#TEST FUNCTIONS

osc1ex = lambda t: np.array([cos(t**2/2), sin(t**2/2)])
def osc1(t, y, k):
    k[0] = -t*y[1]
    k[1] =  t*y[0]

osc2ex = lambda t: np.array([exp(sin(t*t)), exp(cos(t*t))])
def osc2(t, y, k):
    if y[1] > 1e-3:
        k[0] = 2*t*y[0]*log(y[1])
    else:
        k[0] = 2*t*y[0]*log(1e-3)
    if y[0] > 1e-3:
        k[1] = -2*t*y[1]*log(y[0])
    else:
        k[1] = -2*t*y[1]*log(1e-3)

dahlex = lambda t: np.array([exp(-t)])
def dahl(t, y, k):
    k[0] = -y[0]

try:
    brus_sol0 = np.fromfile('brus-sol/sol0')
    brus_sol1 = np.fromfile('brus-sol/sol1')
    brus_t = np.fromfile('brus-sol/t')
except FileNotFoundError:
    print("At least one of the brusselator 'exact' solution files was not found. These files were generated with the 8th order C++ solver at very small time step and copies of them are on my google drive. The exact solution to the brusselator system can't be used without those files, so it is left as None in the test_funks DataFrame.")
    brusex = None
else:
    brus_interp0 = interp1d(brus_t, brus_sol0, kind='linear')
    brus_interp1 = interp1d(brus_t, brus_sol1, kind='linear')
    brusex = lambda t: np.array([brus_interp0(t), brus_interp1(t)])
def brus(t, y, k):
    k[0] = 1 + y[0]*(y[0]*y[1] - 4)
    k[1] = y[0]*(3 - y[0]*y[1])

#-------------------------------------------------------------------------------
#put the functions into a DataFrame

test_funks = pd.DataFrame({
    'name':
        ['dahl',
        'osc1',
        'osc2',
        'brus'],
    'long_name':
        ['Dahlquist Test Equation',
        'Oscillator 1',
        'Oscillator 2',
        'Brusselator'],
    'funk':
        [dahl,
        osc1,
        osc2,
        brus],
    'initial_values':
        [[1],
        [1, 0],
        [1, exp(1)],
        [1.5, 3]],
    'exact_solution':
        [dahlex,
        osc1ex,
        osc2ex,
        brusex]
})

test_funks.set_index('name', inplace=True)
