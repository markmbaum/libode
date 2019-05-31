import sys
from os.path import join
import numpy as np
from numpy import log10
import matplotlib.pyplot as plt
plt.style.use('dark_background')

from test_funks import test_funks
sys.path.append(join('..', '..'))
from pylibode import solvers

#-------------------------------------------------------------------------------
#INPUTS

#name of solver to use (must correspond to a name in the solver DataFrame)
if len(sys.argv) > 1:
    solver = sys.argv[1]
else:
    solver = 'Butcher6'
#name of test system to use
if len(sys.argv) > 2:
    funk = sys.argv[2]
else:
    funk = 'osc2'
#stopping time (tend)
tend = 5
#maximum time step length
dtmax = 0.05
#minimum time step length
dtmin = 0.001
#number of solves
dts = 12
#maximum order line to plot
maxord = 5

#-------------------------------------------------------------------------------
#FUNCTIONS

def conv(solver, funk):

    dt = np.logspace(log10(dtmin), log10(dtmax), dts)[::-1]
    exact = test_funks.at[funk, 'exact_solution']
    err = []
    for i in range(dts):
        #set up the solver
        s = solvers['constructor'][solver](
                test_funks.at[funk, 'funk'],
                test_funks.at[funk, 'initial_values'])
        #solve
        t, sol = s.solve_fixed(tend, dt[i])
        print('finished with dt = %g' % dt[i])
        #compute error
        ex = exact(t)
        err.append(np.sqrt(np.sum((ex - sol)**2)))

    return(dt, err)

def plot_conv(h, err, maxord=6):
    #instantiate plot objs
    fig, ax = plt.subplots(1,1)
    #plot the convergence results, step size and a measure of error
    ax.loglog(h, err, 'w.', label='results')
    #get a dense array between the largest and smallest step size, w/ log spacing
    hh = np.logspace(log10(min(h)), log10(max(h)), 2500)
    #plot ideal convergence lines for orders up to maxord
    for i in range(1, maxord+1):
        #find the vertical shift required to match the largest h (tricky in log)
        b = log10(err[0]/(10**(i*log10(h[0]))))
        #create the line, cutting off portions that extend beyond the results
        y = 10**(i*log10(hh) + b)
        mask = (y >= min(err))
        hh = hh[mask]
        y = y[mask]
        #plot
        ax.loglog(hh, y, ':', label='$\propto (\Delta t)^{%d}$' % i)
        #label
        ax.text(min(hh), min(y), str(i), ha='right', va='bottom')
    ax.set_xlabel('Step Size')
    ax.set_ylabel('$l_2$ error')
    ax.legend()
    return(fig, ax)

#-------------------------------------------------------------------------------
#MAIN

h, err = conv(solver, funk)
fig, ax = plot_conv(h, err, maxord=solvers['order'][solver]+1)
ax.set_title('Convergence of %s' % solvers['long_name'][solver])
plt.show()
