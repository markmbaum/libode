from os.path import join
from numpy import fromfile
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#read the output files
x = fromfile(join('out', 'lorentz_0'))
y = fromfile(join('out', 'lorentz_1'))
z = fromfile(join('out', 'lorentz_2'))
t = fromfile(join('out', 'lorentz_t'))

#plot in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z)
ax.plot([x[0]], [y[0]], [z[0]], 'g.')
ax.plot([x[-1]], [y[-1]], [z[-1]], 'r.')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('lorentz attractor')

#plot cross sections
fig, ax = plt.subplots(1,1)
ax.plot(x,y)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('x-y projection')

fig, ax = plt.subplots(1,1)
ax.plot(y,z)
ax.set_xlabel('y')
ax.set_ylabel('z')
ax.set_title('y-z projection')

fig, ax = plt.subplots(1,1)
ax.plot(x,z)
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_title('x-z projection')

plt.show()
