import numpy as np
import matplotlib.pyplot as plt

#plt.rc('font', family='serif')
#plt.rc('text', usetex=True)

#oscillator 1
#sol1ex = lambda t: np.cos(t**2/2)
#sol2ex = lambda t: np.sin(t**2/2)
#name = 'Osc1'

#oscillator 2
sol1ex = lambda t: np.exp(np.sin(t**2))
sol2ex = lambda t: np.exp(np.cos(t**2))
name = 'Osc2'

t = np.fromfile('../out/%s_t' % name)
sol1 = np.fromfile('../out/%s_0' % name)
sol2 = np.fromfile('../out/%s_1' % name)

fig, axs = plt.subplots(2, 4, figsize=(10,5))
axs = [item for sublist in axs for item in sublist]

tdense = np.linspace(min(t), max(t), 2500)
axs[0].plot(tdense, sol1ex(tdense), 'k', label='$y_1$ exact')
axs[0].plot(t, sol1, '--', label='$y_1$ numerical')
axs[0].set_title('Solutions')
axs[0].set_ylabel('$y_1$')
axs[0].legend()

axs[4].plot(tdense, sol2ex(tdense), 'k', label='$y_2$ exact')
axs[4].plot(t, sol2, '--', label='$y_2$ numerical')
axs[4].set_ylabel('$y_2$')
axs[4].set_xlabel('$t$')
axs[4].legend()

axs[1].semilogy(t, np.abs(sol1 - sol1ex(t)), label='$y_1$ abs err')
axs[5].semilogy(t, np.abs(sol2 - sol2ex(t)), label='$y_2$ abs err')
axs[1].set_title('Absolute Error')
axs[1].set_ylabel('$y_1$ error')
axs[5].set_xlabel('$t$')
axs[5].set_ylabel('$y_2$ error')

axs[2].semilogy(t, np.abs((sol1 - sol1ex(t))/sol1ex(t)), label='$y_1$ rel err')
axs[6].semilogy(t, np.abs((sol2 - sol2ex(t))/sol1ex(t)), label='$y_2$ rel err')
axs[2].set_title('Relative Error')
axs[2].set_ylabel('$y_1$ relative error')
axs[6].set_xlabel('$t$')
axs[6].set_ylabel('$y_2$ relative error')

axs[3].plot(t[:-2], np.diff(t)[:-1])
axs[3].set_ylabel('$\Delta t$')
axs[7].semilogy(t[:-2], np.diff(t)[:-1])
axs[3].set_title('Step Size')

axs[7].set_ylabel('$\log_{10}(\Delta t)$')
axs[7].set_xlabel('$t$')

plt.tight_layout()
plt.show()
