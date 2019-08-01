import numpy as np
import matplotlib.pyplot as plt

#plt.rc('font', family='serif')
#plt.rc('text', usetex=True)

sol1err = np.fromfile('../out/sol1err')
sol2err = np.fromfile('../out/sol2err')
L2err = np.sqrt(sol2err**2 + sol1err**2)
h = np.fromfile('../out/h')
x = np.sort(h)

fig, ax = plt.subplots(1,1)

for i in range(1,10):
    hh = np.logspace(np.log10(min(h)), np.log10(max(h)), 2500)
    b = np.log10(L2err[0]/(10**(i*np.log10(h[0]))))
    y = 10**(i*np.log10(hh) + b)
    mask = (y > min(L2err))
    hh = hh[mask]
    y = y[mask]
    ax.loglog(hh, y, ':', label='$\propto (\Delta t)^{%d}$' % i)
    ax.text(min(hh), min(y), str(i), ha='right', va='bottom')

ax.loglog(h, L2err, 'k.', label='results')
ax.set_xlabel('step size $(\Delta t)$')
ax.set_ylabel('$l_2$ error')
ax.legend()
ax.set_title('Convergence Test')

plt.tight_layout()
plt.show()
