from os.path import join
from numpy import fromfile
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

#load results
q = [fromfile(join('out', 'star_%d' % i)) for i in range(3)]
p = [fromfile(join('out', 'star_%d' % i)) for i in range(3,6)]
t = fromfile(join('out', 'star_t'))
H = fromfile(join('out', 'star_H'))

#make a 3d plot of the trajectory
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(*q, linewidth=0.5)
ax.set_xlabel('$q_1$')
ax.set_ylabel('$q_2$')
ax.set_zlabel('$q_3$')
ax.set_title('Star Trajectory')
fig.tight_layout()

#make cross sections of trajectory
fig, ax = plt.subplots(1,1)
ax.plot(q[0], q[1], linewidth=0.5)
ax.set_xlabel('$q_1$')
ax.set_ylabel('$q_2$')
#ax.axis('equal')
ax.set_title('Projected Star Trajectory')
fig.tight_layout()
fig, ax = plt.subplots(1,1)
ax.plot(q[1], q[2], linewidth=0.5)
ax.set_xlabel('$q_2$')
ax.set_ylabel('$q_3$')
#ax.axis('equal')
ax.set_title('Projected Star Trajectory')
fig.tight_layout()
fig, ax = plt.subplots(1,1)
ax.plot(q[0], q[2], linewidth=0.5)
ax.set_xlabel('$q_1$')
ax.set_ylabel('$q_3$')
#ax.axis('equal')
ax.set_title('Projected Star Trajectory')
fig.tight_layout()

#plot the Hamiltonian over time
fig, ax = plt.subplots(1,1)
ax.plot(t, H)
ax.set_xlabel('$t$')
ax.set_ylabel('$H(q,p)$')
ax.set_title('Hamiltonian (conservation test)')
fig.tight_layout()

for i in range(len(q)):
    q[i] = q[i][::3]

fig, (axa, axb, axc) = plt.subplots(3,1)
axa.set_xlabel('$q_1$')
axa.set_ylabel('$q_2$')
axa.set_xlim(min(q[0]), max(q[0]))
axa.set_ylim(min(q[1]), max(q[1]))
axb.set_xlabel('$q_2$')
axb.set_ylabel('$q_3$')
axb.set_xlim(min(q[1]), max(q[1]))
axb.set_ylim(min(q[2]), max(q[2]))
axc.set_xlabel('$q_1$')
axc.set_ylabel('$q_3$')
axc.set_xlim(min(q[0]), max(q[0]))
axc.set_ylim(min(q[2]), max(q[2]))
la, = axa.plot([], [])
lb, = axb.plot([], [])
lc, = axc.plot([], [])

def init():
    return(la, lb, lc)

def animate(i):
    la.set_xdata(q[0][:i])
    la.set_ydata(q[1][:i])
    lb.set_xdata(q[1][:i])
    lb.set_ydata(q[2][:i])
    lc.set_xdata(q[0][:i])
    lc.set_ydata(q[2][:i])
    return(la, lb, lc)

fig.tight_layout()
ani = animation.FuncAnimation(fig, animate, init_func=init, blit=True, interval=0.01, frames=len(q[0]))

plt.show()
