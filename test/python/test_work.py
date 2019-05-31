import sys
from os.path import join
import numpy as np
from time import time
import matplotlib.pyplot as plt
plt.style.use('dark_background')

from test_funks import test_funks
sys.path.append(join('..', '..'))
from pylibode import solvers

#-------------------------------------------------------------------------------
#INPUTS

#name of solver to use (must correspond to a name in the solver DataFrame)
solver_names = ['Euler', 'Ssp3', 'RK4', 'DoPri54', 'DoPri87']
#name of test system to use
if len(sys.argv) > 2:
    funk = sys.argv[2]
else:
    funk = 'osc2'
#time step boundaries to use for each (in log10 units)
dt_bounds = [(-3,-5), (-1,-3), (-1,-3), (-1,-3), (-0.5,-1.5)]
#number of time step intervals
ndt = 6
#stopping time (tend)
tend = 10

#-------------------------------------------------------------------------------
#FUNCTIONS

def measure_work(solver, dt1, dt2):

    ex = test_funks['exact_solution'][funk]
    dts = np.logspace(dt1, dt2, ndt)
    ts = np.empty((ndt,))
    err = np.empty((ndt,))
    neval = np.empty((ndt,), dtype=int)
    for i in range(ndt):
        s = solvers['constructor'][solver](
                test_funks['funk'][funk],
                test_funks['initial_values'][funk])
        t1 = time()
        t, sol = s.solve_fixed(tend, dts[i])
        t2 = time()
        ts[i] = t2 - t1
        err[i] = np.sqrt(np.sum((sol - ex(t))**2))
        neval[i] = s.neval
    return(dts, ts, err, neval)

#-------------------------------------------------------------------------------
#MAIN

fig, (ax1, ax2) = plt.subplots(1,2)
for dtb, solver_name in zip(dt_bounds, solver_names):
    print(solver_name)
    dts, ts, err, neval = measure_work(solver_name, *dtb)
    ax1.loglog(err, neval, '.:', label=solver_name)
    ax2.loglog(err, ts, '.:', label=solver_name)

ax1.legend()
ax1.set_xlabel("$L_2$ error")
ax1.set_ylabel('Function Evaluations')

ax2.legend()
ax2.set_xlabel("$L_2$ error")
ax2.set_ylabel('Execution Time ($s$)')

fig.tight_layout()
plt.show()
