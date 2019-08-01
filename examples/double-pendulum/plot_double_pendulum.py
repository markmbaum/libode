from os.path import join
import numpy as np
from numpy import pi, nan, full, array, fromfile
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import LineCollection

#load results
x1 = fromfile(join('out','x1'))
y1 = fromfile(join('out','y1'))
x2 = fromfile(join('out','x2'))
y2 = fromfile(join('out','y2'))
t1 = fromfile(join('out','double_pendulum_0'))
t2 = fromfile(join('out','double_pendulum_2'))
t = fromfile(join('out','double_pendulum_t'))

L1 = (x1[0]**2 + y1[0]**2)**(0.5)
L2 = ((x1[0] - x2[0])**2 + (y1[0] - y2[0])**2)**(0.5)
L = (L1 + L2)*1.1

fig = plt.figure()
axa = plt.subplot2grid((2,4), (0,0), colspan=2, rowspan=2)
axa.set_xlim(-L, L)
axa.set_ylim(-L, L)
axa.set_aspect('equal')
axa.set_xlabel('$x$')
axa.set_ylabel('$y$')
axa.set_title('Double Pendulum')

line1, = axa.plot([0,x1[0]], [0,y1[0]], ':', color='gray', zorder=-1, alpha=0.5)
dot1, = axa.plot(x1[0], y1[0], '.', color='gray', markersize=15, zorder=1)
trace1, = axa.plot(x1[:0], y1[:0], color='gray', linewidth=0.75, alpha=0.5)

line2, = axa.plot([x1[0], x2[0]], [y1[0],y2[0]], 'k:', zorder=-1, alpha=0.5)
dot2, = axa.plot(x2[0], y2[0], 'k.', markersize=15, zorder=1)
verts = list(zip(x2, y2))
segs = [[verts[i],verts[i+1]] for i in range(len(verts)-1)]
trace2 = LineCollection(segs)
colors = full((len(x1),), 'none')
trace2.set_color(colors)
trace2.set_alpha(0.8)
axa.add_collection(trace2)

axb = plt.subplot2grid((2,4), (0,2), colspan=2, rowspan=1)
axb.set_xlim(-0.01*t.max(), 1.01*t.max())
axb.set_ylim(min(t1)*1.01, max(t1)*1.01)
axb.set_ylabel(r'$\theta_1$')
th1, = axb.plot([], [], color='gray')

axc = plt.subplot2grid((2,4), (1,2), colspan=2, rowspan=1)
axc.set_xlim(-0.01*t.max(), 1.01*t.max())
axc.set_ylim(min(t2)*1.01, max(t2)*1.01)
axc.set_ylabel(r'$\theta_2$')
axc.set_xlabel('time')
th2, = axc.plot([], [], 'k')

def init():
    line1.set_ydata(array(line1.get_ydata())*nan)
    dot1.set_ydata(array(dot1.get_ydata())*nan)
    trace1.set_ydata(array(trace1.get_ydata())*nan)
    line2.set_ydata(array(line2.get_ydata())*nan)
    dot2.set_ydata(array(dot2.get_ydata())*nan)
    return(line1, dot1, trace1, line2, dot2, trace2, th1, th2)

def animate(i):
    line1.set_xdata([0,x1[i]])
    line1.set_ydata([0,y1[i]])
    dot1.set_xdata([x1[i]])
    dot1.set_ydata([y1[i]])
    trace1.set_xdata(x1[:i])
    trace1.set_ydata(y1[:i])

    line2.set_xdata([x1[i],x2[i]])
    line2.set_ydata([y1[i],y2[i]])
    dot2.set_xdata([x2[i]])
    dot2.set_ydata([y2[i]])

    j = i - 400
    if j < 0: j = 0
    colors[:] = 'none'
    colors[j:i] = 'k'
    trace2.set_color(colors)

    th1.set_xdata(t[:i])
    th1.set_ydata(t1[:i])
    th2.set_xdata(t[:i])
    th2.set_ydata(t2[:i])

    return(line1, dot1, trace1, line2, dot2, trace2, th1, th2)

ani = animation.FuncAnimation(fig, animate, init_func=init, blit=True, interval=0.01, frames=len(x1))

fig.tight_layout()
plt.show()
