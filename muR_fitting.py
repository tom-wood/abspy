import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import chebval, chebfit
from numpy.polynomial.polynomial import polyval

muR_data = np.genfromtxt('C:\\Users\\vgx18551\\Documents\\Data\\I11_Feb17_TMnitriding\\muR_intensities.txt', 
                         skip_header=1)

"""
#initial plot of data
fig = plt.figure()
ax = fig.add_subplot(111)
labels=[str(n * 5) for n in range(19)]
for i in range(muR_data.shape[1] - 1):
    ax.plot(muR_data[:, 0], muR_data[:, i + 1], label=labels[i])
ax.legend(loc=0)
fig.tight_layout()
"""

#Take values for muR=18.1
i18_1 = np.searchsorted(muR_data[:, 0], 18.1)

poly_rank = 5
x_data = np.array([n * 10 for n in range(19)]) #in 2 theta in degrees
x_data = x_data * (np.pi / 180.)
x_fit = np.linspace(x_data[0], x_data[-1], 1000)
pvals = np.zeros((muR_data.shape[0] - 1, poly_rank + 1))
y_data = np.zeros((muR_data.shape[0] - 1, muR_data.shape[1] - 1))
yfits = np.zeros((muR_data.shape[0] - 1, 1000))
for i in range(muR_data.shape[0] - 1):
    yvals = muR_data[i + 1, 1:]
    yvals = yvals / yvals.max()
    p = np.polyfit(x_data, yvals, poly_rank)
    pvals[i] = p
    y_fit = polyval(x_fit, p[::-1])
    yfits[i] = y_fit
    
##Now fit the pvalues
pfit_order = 16
pfits = np.zeros((poly_rank + 1, pfit_order + 1))
pfit_x = np.linspace(muR_data[1, 0], muR_data[-1, 0], 1000)
pfit_y = np.zeros((poly_rank + 1, 1000))
muR_mean = muR_data[:, 0].mean()
for i in range(poly_rank + 1):
    #pfit = chebfit(muR_data[:, 0], pvals[:, i], pfit_order - 1)
    pfit = np.polyfit(muR_data[1:, 0], pvals[:, i], pfit_order)
    pfits[i] = pfit
    pfit_y[i] = polyval(pfit_x, pfit[::-1])
#pfit_5 = np.polyfit(muR_data[1:, 0], pvals[:, 4], 12)
#pfits[4, 5:] = pfit_5
#pfit_y[4] = polyval(pfit_x, pfit_5[::-1])

#Now try plotting these values
#for i in range(poly_rank + 1):
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.plot(muR_data[1:, 0], pvals[:, i], 'go')
#    ax.plot(pfit_x, pfit_y[i], 'g')
#    ax.set_xlabel(r'$\mu$R')
#    ax.set_ylabel('p values')
#    fig.tight_layout()

#Now reverse the fits as a check...
rev_fits_p = np.zeros(pvals.shape).T
for i in range(pfits.shape[0]):
    rev_fits_p[i] = polyval(muR_data[1:, 0], pfits[i][::-1])

rev_fits_y = np.zeros((yfits.shape[0], x_data.shape[0]))
for i in range(yfits.shape[0]):
    rev_fits_y[i] = polyval(x_data, rev_fits_p[:, i][::-1])

fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(muR_data.shape[0] - 1):
    if i % 10 == 0:#i == len(range(muR_data.shape[0] - 1)) - 1:
        yvals = muR_data[i + 1, 1:]
        yvals = yvals / yvals.max()
        ax.plot(x_data, yvals, 'o')
        ax.plot(x_data, rev_fits_y[i])
ax.set_xlabel(r'$\theta$')
ax.set_ylabel('Intensity')
fig.tight_layout()

#Plot p values from muR >= 5:
#i5 = np.searchsorted(muR_data[:, 0], 5)
#pfit_order = 10
#pfits = np.zeros((5, pfit_order))
#pfit_x = np.linspace(muR_data[i5, 0], muR_data[-1, 0], 1000)
#pfit_y = np.zeros((5, 1000))
##muR_mean = muR_data[:, 0].mean()
#for i in range(5):
#    #pfit = chebfit(muR_data[:, 0], pvals[:, i], pfit_order - 1)
#    pfit = np.polyfit(muR_data[i5:, 0], pvals[i5 - 1:, i], pfit_order - 1)
#    pfits[i] = pfit
#    pfit_y[i] = polyval(pfit_x, pfit[::-1])
#for i in range(5):
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.plot(muR_data[i5:, 0], pvals[i5 - 1:, i], 'go')
#    ax.plot(pfit_x, pfit_y[i], 'g')
#    ax.set_xlabel(r'$\mu$R')
#    ax.set_ylabel('p values')
#    fig.tight_layout()

#And plot the fits
#fig = plt.figure()
#ax = fig.add_subplot(111)
#for i in range(muR_data.shape[0] - 1):
#    if i % 10 == 0:#i == len(range(muR_data.shape[0] - 1)) - 1:
#        yvals = muR_data[i + 1, 1:]
#        yvals = yvals / yvals.max()
#        ax.plot(x_data, yvals, 'o')
#        ax.plot(x_fit, yfits[i])
#ax.set_xlabel(r'$\theta$')
#ax.set_ylabel('Intensity')
#fig.tight_layout()


