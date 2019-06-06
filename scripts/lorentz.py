import sys
sys.path.append('..')

from pylibode import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#define the Lorentz system
def ode_funk(t, y, k):
    k[0] = sigma*(y[1] - y[0])      # dy_0 / dt
    k[1] = y[0]*(rho - y[2]) - y[1] # dy_1 / dt
    k[2] = y[0]*y[1] - beta*y[2]    # dy_2 / dt

#define the system's parameters
sigma = 10
beta = 8/3
rho = 28

#create a solver and solve with a fixed time step
solver = solvers['constructor']['DoPri87'](ode_funk, [1.]*3)
t, sol = solver.solve_fixed(20, 0.005, snaps='all')
x, y, z = sol[:,0], sol[:,1], sol[:,2]

#plot in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
