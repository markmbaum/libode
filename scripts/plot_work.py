import numpy as np
import matplotlib.pyplot as plt

#plt.rc('font', family='serif')
#plt.rc('text', usetex=True)

def plot_work(name):
    sol1e = np.fromfile('../out/sol1err_' + name)
    sol2e = np.fromfile('../out/sol2err_' + name)
    f = np.fromfile('../out/neval_' + name)
    ax1.loglog(sol1e, f, label=name)
    ax2.loglog(sol2e, f, label=name)
    ax3.loglog(np.sqrt(sol2e**2 + sol1e**2), f, label=name)


fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,5))
fig2, ax3 = plt.subplots(1, 1, figsize=(8,5))

plot_work('Euler')
plot_work('Trapz')
plot_work('Ssp3')
plot_work('RK4')
plot_work('DoPri54')
plot_work('Vern65')
plot_work('Vern76')
plot_work('DoPri87')
plot_work('Vern98')

ax1.set_title('$y_1$')
ax1.set_xlabel('Error')
ax1.set_ylabel('Function Evaluations')
ax1.invert_xaxis()
ax1.legend()

ax2.set_title('$y_2$')
ax2.set_xlabel('Error')
ax2.invert_xaxis()
ax2.legend()

ax3.set_title('Accuracy vs Function Calls for Different Methods')
ax3.set_xlabel('Precision ($l_2$ Error)')
ax3.set_ylabel('Function Evaluations')
ax3.invert_xaxis()
#box = ax3.get_position()
#ax3.set_position([box.x0, box.y0, box.width * 0.85, box.height])
ax3.legend()#loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()
