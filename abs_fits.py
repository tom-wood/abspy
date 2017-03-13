import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import pandas as pd
from numba import jit

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

def fit_to_cosh_sinh(y, fy, guesses, fit_array=True):
    """Return best fit for y vs fy with sum of sinh and cosh"""
    if fit_array:
        y_fit = np.linspace(-1, 1, 1000)
        a, b, c, d, e = guesses
        fy_guess = a * np.sinh(b * y_fit) + c * np.cosh(d * y_fit) + e
    def residuals(guesses, y, fy):
        a, b, c, d, e = guesses
        fy_calc = a * np.sinh(b * y) + c * np.cosh(d * y) + e
        err = np.abs(fy - fy_calc)
        return err
    p = list(leastsq(residuals, guesses, args=(y, fy)))[0]
    if fit_array:
        a, b, c, d, e = p
        fy_fit = a * np.sinh(b * y_fit) + c * np.cosh(d * y_fit) + e
        return p, y_fit, fy_fit, fy_guess
    return p

########this one doesn't work very well##########
@jit(nopython=True)
def res_hyp(guesses, y, fy):
    a, b, c, d, e = guesses
    fy_calc = (a * np.sinh(b * y) + c * np.cosh(d * y)) * e * \
              np.sqrt((1 - y**2))
    err = np.sum(np.abs(fy - fy_calc))
    return err

########this one doesn't work very well either##########
@jit(nopython=True)
def res_exp(guesses, y, fy):
    a, b, c, d, e = guesses
    fy_calc = (a * np.sinh(b * y) + c * np.cosh(d * y)) * e * \
              np.sqrt((1 - y**2))
    err = np.sum(np.abs(fy - fy_calc))
    return err

@jit(nopython=True)
def res_exp_con(guesses, y, fy):
    a, b, c, d, e, f = guesses
    fy_calc = (a * np.sinh(b * y) + c * np.cosh(d * y) + f) * e * \
              np.sqrt((1 - y**2))
    err = np.sum(np.abs(fy - fy_calc))
    return err

#Now try putting into polar coordinates ####doesn't help
def get_polar(y, fy):
    """Puts y, fy into polar coordinates"""
    dists = np.sqrt(y**2 + fy**2)
    angles = np.arctan2(fy, y)
    return dists, angles
