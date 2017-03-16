import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pandas as pd
from numba import jit

def fig_setup():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    return fig, ax

muvals = np.arange(0.1, 20.05, 0.1) 
thetavals = np.arange(0, 90.5, 1) * np.pi / 180.
num_yvals = 400
yvals = np.linspace(-1., 1., num_yvals)
ydiff = yvals[1] - yvals[0]
fy = np.zeros((len(thetavals), len(muvals), len(yvals)))

fnames = ['abs_1mmdiam_capillary_theta' + str(r) + '.csv' for r in 
          range(91)]
for i, fn in enumerate(fnames):
    data = pd.read_csv(fn, sep=',')
    fy[i] = data.values[:, 1:].T
print 'ok'

#mu_i = np.searchsorted(muvals, 18)
#mu_fits = np.zeros((fy.shape[0], poly_deg + 1))
#for th_i in range(fy.shape[0]):
#    p = np.polyfit(yvals[1:-1], fy[th_i, mu_i, 1:-1] / \
#                                np.sqrt(1 - yvals[1:-1]**2), 20)
#    mu_fits[th_i] = p

def get_offsets(yvals, fy):
    offsets = np.zeros(fy.shape[:2])
    for th_i in range(fy.shape[0]):
        for mu_i in range(fy.shape[1]):
            cum_fy = np.cumsum(fy[th_i, mu_i])
            #centre_i = np.searchsorted(cum_fy, cum_fy[-1] / 2.)
            #offsets[th_i, mu_i] = yvals[centre_i]
            offsets[th_i, mu_i] = np.interp(cum_fy[-1] / 2., cum_fy, yvals)
    return offsets

offsets = get_offsets(yvals, fy)

#now fit to muval of 3.2
#mu_i = np.searchsorted(muvals, 3.2)
#mu_i = -1
poly_deg = 6
off_ps = np.zeros((offsets.shape[1], poly_deg + 1))
for mu_i in range(len(muvals)):
    off_ps[mu_i] = np.polyfit(thetavals[1:], offsets[1:, mu_i], poly_deg)

#th_fit = np.linspace(0, thetavals[-1], 1000)
#off_fit = np.polyval(off_p, th_fit)
#
#fig, ax = fig_setup()
#ax.plot(thetavals, offsets[:, mu_i], 'o')
#ax.plot(th_fit, off_fit)

#next fit the polynomial factors themselves to mu

p_deg = 18
p_fits = np.zeros((poly_deg + 1, p_deg + 1))
for i in range(poly_deg + 1):
    p_fits[i] = np.polyfit(muvals, off_ps[:, i], p_deg)

#p_fits_x = np.linspace(muvals[0], muvals[-1], 1000)
#for i in range(poly_deg + 1):
#    fig, ax = fig_setup()
#    ax.plot(muvals, off_ps[:, i], 'o')
#    ax.plot(p_fits_x, np.polyval(p_fits[i], p_fits_x))

p_abs = 'If_Prm_Eqn_Rpt(murc, murv, min 0.1, max 20)\n'
for i1 in range(poly_deg + 1):
    p_str = 'prm #m_unique p_off' + str(i1) + ' = '
    for i2 in range(p_deg + 1):
        p_str = p_str + '(' + str(p_fits[i1][i2]) + \
                ') * CeV(murc, murv)^%d + ' % (p_deg - i2)
    p_str = p_str[:-2]
    p_str += ';\n'
    p_abs += p_str

p_abs2 = 'th2_offset = ('
for i1 in range(poly_deg + 1):
    p_abs2 = p_abs2 + 'p_off' + str(i1) + '* Th^%d + ' % (poly_deg - i1)

p_abs2 = p_abs2[:-3]
p_abs2 += ') * (0.5 / Rs) Rad;' 

tstr = 'macro Cylinder_Peak_Position {\n' 
tstr += '#m_argu murc\n'
tstr += p_abs + p_abs2 + '}'

write = True
if write:
    with open('macro_th2offset3pt1.txt', 'w') as f:
        f.write(tstr)

