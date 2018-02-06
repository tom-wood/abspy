import numpy as np
import matplotlib.pyplot as plt
from numba import jit

#start off for mu = 1, r = 0.5, theta = 0 (because, well, why not?)
r = 0.35
save_fname = 'abs_pt7mmdiam_capillary_theta'
muvals = np.arange(0.1, 20.05, 0.1) #these are actually mu*r values 
thetavals = np.arange(0, 90.5, 1) * np.pi / 180.
num_yvals = 400
write_vals = True

@jit(nopython=True)
def get_vals(muvals, thetavals, r, num_yvals):
    yvals = np.linspace(-1., 1., num_yvals)
    ydiff = yvals[1] - yvals[0]
    fy = np.zeros((len(muvals), len(thetavals), len(yvals)))
    for im, mu in enumerate(muvals):
        for it, theta in enumerate(thetavals):
            for i1, y in enumerate(yvals):
                zmax = np.sqrt(1 - y**2)
                num_zvals = int(2 * zmax / ydiff) + 1
                zvals = np.linspace(-zmax, zmax, num_zvals)
                fybit = np.zeros(len(zvals) - 1)
                for i2, z in enumerate(zvals[1:]):
                    x1 = r * (z * np.cos(2 * theta) - y * np.sin(2 * theta)\
                              + np.sqrt(1 - (y * np.cos(2 * theta) + \
                                             z * np.sin(2 * theta))**2))
                    x2 = r * (-z + np.sqrt(1 - y**2))
                    fybit[i2] = np.exp(-mu * (x1 + x2)) * (z - zvals[i2])
                fy[im, it, i1] = 0.5 * np.sum(fybit)
        print mu
    return fy

fy = get_vals(muvals, thetavals, r, num_yvals)

if write_vals:
    fnames = [save_fname + str(th) + '.csv' for th in 
              range(91)]
    muvals_a = np.concatenate((np.array([0]), muvals))
    yvals = np.linspace(-1, 1., num_yvals)
    yvals = yvals.reshape(yvals.shape[0], 1)
    for i in range(91):
        values = np.column_stack((yvals, fy[:, i, :].T))
        values = np.row_stack((muvals_a, values))
        np.savetxt(fnames[i], values, delimiter=',')
