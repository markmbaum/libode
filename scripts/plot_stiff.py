import numpy as np
import matplotlib.pyplot as plt

#plt.rc('font', family='serif')
#plt.rc('text', usetex=True)

name = "Brus"

t = np.fromfile('../out/%s_t' % name)
sol1 = np.fromfile('../out/%s_0' % name)
sol2 = np.fromfile('../out/%s_1' % name)

fig, axs = plt.subplots(3, 1, figsize=(10,5))

axs[0].plot(t, sol1)
axs[0].set_ylabel('$y_1$')

axs[1].plot(t, sol2)
axs[1].set_ylabel('$y_2$')

axs[2].plot(t[:-2], np.diff(t)[:-1])
axs[2].set_ylabel('$\Delta t$')
axs[2].set_xlabel('$t$')

plt.tight_layout()
plt.show()
