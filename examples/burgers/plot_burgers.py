from os import listdir
from os.path import join
from numpy import fromfile
import matplotlib.pyplot as plt

#read the output files
fns = [fn for fn in listdir('out') if ('_t' not in fn and fn != 'x')]
fns = sorted(fns, key=lambda fn: int(fn.split('_')[-1]))
u = [fromfile(join('out', fn)) for fn in fns]

t = fromfile(join('out', 'burgers_snap_t'))
x = fromfile(join('out', 'x'))
print('%d spatial nodes' % len(x))

fig, ax = plt.subplots(1,1)
ax.plot(t, [sum(u[i])*(x[1] - x[0]) for i in range(len(u))], 'ko')
ax.set_xlabel('$t$')
ax.set_ylabel('Total $u$')
ax.set_title('Total $u$ (conservation test)')
fig.tight_layout()

fig, ax = plt.subplots(1,1)
N = len(u)
for i in range(len(u)):
    ax.plot(x, u[i], 'k', alpha=(0.1 + 0.9*(i/N)), label='$t=%g$' % t[i])
#ax.legend()
ax.set_xlabel('$x$')
ax.set_ylabel('$u$')
ax.set_title("Viscid Burger's Equation\n$u_t = -u \cdot u_x + u_{xx}/1000$")
fig.tight_layout()

plt.show()
