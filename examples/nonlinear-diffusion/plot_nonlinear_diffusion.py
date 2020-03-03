from numpy import *
from os import listdir
from os.path import join
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

#read the output files
fns = [fn for fn in listdir('out') if ('_t' not in fn and fn != 'x')]
fns = sorted(fns, key=lambda fn: int(fn.split('_')[-1]))
u = [fromfile(join('out', fn)) for fn in fns]

t = fromfile(join('out', 'nonlinear_diffusion_snap_t'))
x = fromfile(join('out', 'x'))
print('%d spatial nodes' % len(x))

fig, ax = plt.subplots(1,1)
ax.plot(t, [sum(u[i])*(x[1] - x[0]) for i in range(len(u))], 'ko')
ax.set_xlabel('$t$')
ax.set_ylabel('Total $y$')
ax.set_title('Total $y$ (conservation test)')
fig.tight_layout()

fig, ax = plt.subplots(1,1)
y = array([max(abs(diff(i))) for i in u])
y /= y.max()
ax.plot(t, y, 'k.')
tt = linspace(t.min(), t.max(), 1000)
yy = CubicSpline(t, y)(tt)
ax.plot(tt, yy, 'k', zorder=-1)
ax.set_ylabel('Maximum Slope Magnitude (Normalized)')
ax.set_xlabel('$t$')

fig, ax = plt.subplots(1,1)
N = len(u)
for i in range(len(u)):
    ax.plot(x, u[i], 'k', alpha=(0.1 + 0.9*(i/N)), label='$t=%g$' % t[i])
#ax.legend()
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_title("Nonlinear Diffusion Equation\n$y_t =  y^{3} \cdot y_{xx}$")
fig.tight_layout()

plt.show()
