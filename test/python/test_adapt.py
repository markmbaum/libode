import sys
from os.path import join
import numpy as np
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
    solver = 'DoPri87'
#name of test system to use
if len(sys.argv) > 2:
    funk = sys.argv[2]
else:
    funk = 'osc2'
#stopping time (tend)
tend = 10
#error tolerances
abstol = 1e-6
reltol = 1e-6

#-------------------------------------------------------------------------------
#MAIN

#set up the solver
s = solvers['constructor'][solver](
        test_funks['funk'][funk],
        test_funks['initial_values'][funk])
s.abstol = abstol
s.reltol = reltol

#solve
t, sol = s.solve_adaptive(tend, tend/100, snaps='all')
print('%d steps, %d function evaluations, %d rejected steps'
        % (s.nstep, s.neval, s.nrej))

#plot
fig, (ax1, ax2, ax3) = plt.subplots(3,1)
tt = np.linspace(0, tend, 2500)
exact = test_funks.at[funk, 'exact_solution']

# exact solution and computed solution in upper panel
for i in range(sol.shape[1]):
    ax1.plot(tt, [exact(j)[i] for j in tt], 'w', linewidth=0.75)
if s.nstep > 150:
    for i in range(sol.shape[1]):
        ax1.plot(t, sol[:,i])
else:
    for i in range(sol.shape[1]):
        ax1.plot(t, sol[:,i], '.')
ax1.set_ylabel('Solutions')

#L2 error in the lower panel
ex = [exact(t[i]) for i in range(sol.shape[0])]
L2 = [np.sqrt(np.sum((ex[i] - sol[i,:])**2)) for i in range(sol.shape[0])]
if s.nstep > 150:
    ax2.semilogy(t, L2, 'C3')
else:
    ax2.semilogy(t, L2, 'C3.')
ax2.set_ylabel('$L_2$ Relative Error')

#step size
if s.nstep > 150:
    ax3.semilogy(t[:-2], np.diff(t)[:-1], 'w')
else:
    ax3.semilogy(t[:-2], np.diff(t)[:-1], 'w.')
ax3.set_ylabel("$\Delta t$")
ax3.set_xlabel('time')

ax1.set_title(solvers.at[solver, 'long_name'])
fig.tight_layout()
plt.show()
