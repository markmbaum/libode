import numpy as np
import matplotlib.pyplot as plt

flatten = lambda L: [item for sublist in L for item in sublist]

frgb = lambda theta: 0.45*(1 + np.cos(theta))

def rgb(theta):
    r = frgb(theta)
    g = frgb(theta - 2*np.pi/3)
    b = frgb(theta + 2*np.pi/3)
    return(r,g,b)

tsnap = np.fromfile('../out/snap_times')
nsnap = len(tsnap)
snaps = [np.fromfile('../out/snap_%d' % i) for i in range(nsnap)]
n = int(len(snaps[0])/3)

fig, axs = plt.subplots(3,3)
axs = flatten(axs)
for i,ax in enumerate(axs):
    snap = snaps[i]
    x, y, theta = snap[:n], snap[n:2*n], snap[2*n:]
    ax.scatter(x, y, s=2, c=np.array([rgb(z) for z in theta]))
    ax.axis('equal')
    ax.set_title('t = %g' % tsnap[i])

plt.tight_layout()
plt.show()
