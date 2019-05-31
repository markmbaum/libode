import os
import imageio
import numpy as np
import matplotlib
import multiprocessing
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

#plot config
plt.style.use('dark_background')
plt.rc('font', family='monospace')
matplotlib.rcParams['figure.dpi'] = 150
matplotlib.rcParams['savefig.dpi'] = 150

#---------------------------------------------------------------------------
#INPUT

dirsnaps = os.path.join('..', 'out')

fnmov = os.path.join('..', 'swarm_movie.mp4')

fps = 60

cpus = 4

#---------------------------------------------------------------------------
#FUNKS

frgb = lambda theta: 0.45*(1 + np.cos(theta))

def rgb(theta):
    r = frgb(theta)
    g = frgb(theta - 2*np.pi/3)
    b = frgb(theta + 2*np.pi/3)
    return(r,g,b)

def crop(im, mult=16):
    n, m = im.shape[:2]
    dn = n - 16*(n//16)
    dnl = dn//2
    dnr = dn - dnl
    dm = m - 16*(m//16)
    dml = dm//2
    dmr = dm - dml
    im = im[dnl:-dnr, dml:-dmr, ...]
    return(im)

def make_image(snap, n, t, xlim, ylim):
    fig, ax = plt.subplots(1,1)
    ax.set_aspect('equal', 'box')
    ax.set_xlim(min([xlim[0], ylim[0]]), max([xlim[1], ylim[1]]))
    ax.set_ylim(min([xlim[0], ylim[0]]), max([xlim[1], ylim[1]]))
    ax.grid(False)
    x, y = snap[:n], snap[n:2*n]
    ax.scatter(x, y, s=1, c=np.array([rgb(z) for z in snap[2*n:]]))
    ax.set_title('Time = %7.2f' % t)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    fig.canvas.draw()
    im = np.array(fig.canvas.renderer._renderer)
    plt.close(fig)
    return(im)


#---------------------------------------------------------------------------
#MAIN

if __name__ == '__main__':

    #load results
    tsnap = np.fromfile('../out/snap_times')
    nsnap = len(tsnap)
    snaps = [np.fromfile('../out/snap_%d' % i).astype('float32') for i in range(nsnap)]
    n = int(len(snaps[0])/3)
    minx = min([np.min(i[:n]) for i in snaps])
    maxx = max([np.max(i[:n]) for i in snaps])
    xlim = (minx, maxx)
    miny = min([np.min(i[n:2*n]) for i in snaps])
    maxy = max([np.max(i[n:2*n]) for i in snaps])
    ylim = (miny, maxy)

    #make plots and dump them into a movie on the way
    #P = multiprocessing.Pool(cpus)
    kw = dict(mode='I', fps=fps)
    with imageio.get_writer(fnmov, **kw) as writer:
        i = 0
        while i < nsnap:
            if i + cpus >= nsnap:
                R = range(i,nsnap)
            else:
                R = range(i, i+cpus)

            args = [(snaps[j], n, tsnap[j], xlim, ylim) for j in R]
            #frames = P.starmap(make_image, args)
            frames = [make_image(*a) for a in args]

            for r,im in zip(R, frames):
                print('frame %d' % r)
                writer.append_data(im)

            i += cpus
