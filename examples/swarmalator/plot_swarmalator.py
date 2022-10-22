from os import listdir
from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.style.use('dark_background')
#plt.rc('font', family='monospace')

frgb = lambda theta: 0.45*(1 + np.cos(theta))

def rgb(theta):
    r = frgb(theta)
    g = frgb(theta - 2*np.pi/3)
    b = frgb(theta + 2*np.pi/3)
    return(r,g,b,1)

#load results
t = np.fromfile(join('out', 'swarmalator_t'))
fns = [join('out',fn) for fn in listdir('out') if fn[-1] != 't']
fns = sorted(fns, key=lambda fn: int(fn.split('_')[-1]))
out = [np.fromfile(fn).astype('float16') for fn in fns]
nag = len(out)//3
x, y, theta = out[:nag], out[nag:nag*2], out[nag*2:]
x, y, theta = np.vstack(x), np.vstack(y), np.vstack(theta)
nout = x.shape[1]

def setup_plot(title='Oscillators That Sync and Swarm'):
    fig, ax = plt.subplots(1,1)
    ax.set_xlim(x.min()*1.05, x.max()*1.05)
    ax.set_ylim(y.min()*1.05, y.max()*1.05)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title)
    return fig, [ax.plot(x[i,0], y[i,0], '.', color=rgb(theta[i,0]), markersize=10)[0] for i in range(nag)]

def init():
    return(L)

def animate(i):
    for j,l in enumerate(L):
        l.set_xdata(x[j,i])
        l.set_ydata(y[j,i])
        l.set_color(rgb(theta[j,i]))
    return(L)

fig, L = setup_plot()
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(L), interval=0.01, blit=True)
plt.show()

#x, y, theta = x[:,::2], y[:,::2], theta[:,::2]
fig, L = setup_plot('')
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(L), interval=1, blit=True)
ani.save('swarmalator.gif', dpi=150)