import numpy as np
import matplotlib.pyplot as plt

#Dahlquist test
#sol1ex = lambda t: np.exp(-t)
#sol2ex = lambda t: np.exp(-2*t)
#oscillator 1
sol1ex = lambda t: np.cos(t**2/2)
sol2ex = lambda t: np.sin(t**2/2)
#oscillator 2
#sol1ex = lambda t: np.exp(np.sin(t**2))
#sol2ex = lambda t: np.exp(np.cos(t**2))

name = 'Osc1'
t = np.fromfile('../out/%s_snap_t' % name)
nsnap = len(t)
sol1 = np.zeros((nsnap,))
sol2 = sol1.copy()
for i in range(nsnap):
    s = np.fromfile('../out/%s_snap_%d' % (name,i))
    sol1[i] = s[0]
    sol2[i] = s[1]

fig, axs = plt.subplots(2, 3, figsize=(10,5))
axs = [item for sublist in axs for item in sublist]

tdense = np.linspace(min(t), max(t), 2500)
axs[0].plot(tdense, sol1ex(tdense), 'k', linewidth=0.5, label='$y_1$ exact')
axs[0].plot(t, sol1, 'C0.', label='$y_1$ numerical')
axs[0].set_title('Solutions')
axs[0].set_ylabel('$y_1$')
axs[0].legend()

axs[3].plot(tdense, sol2ex(tdense), 'k', linewidth=0.5, label='$y_2$ exact')
axs[3].plot(t, sol2, 'C1.', label='$y_2$ numerical')
axs[3].set_ylabel('$y_2$')
axs[3].legend()

axs[1].semilogy(t, np.abs(sol1 - sol1ex(t)), 'C0.', label='$y_1$ abs err')
axs[4].semilogy(t, np.abs(sol2 - sol2ex(t)), 'C1.', label='$y_2$ abs err')
axs[1].set_title('Absolute Error')

axs[2].semilogy(t, np.abs((sol1 - sol1ex(t))/sol1ex(t)), 'C0.', label='$y_1$ rel err')
axs[5].semilogy(t, np.abs((sol2 - sol2ex(t))/sol1ex(t)), 'C1.', label='$y_2$ rel err')
axs[2].set_title('Relative Error')

axs[3].set_xlabel('t')
axs[4].set_xlabel('t')
axs[5].set_xlabel('t')

plt.tight_layout()
plt.show()
