import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import pandas as pd

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

