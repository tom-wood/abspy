import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
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
    p = list(opt.leastsq(residuals, guesses, args=(y, fy)))[0]
    if fit_array:
        a, b, c, d, e = p
        fy_fit = a * np.sinh(b * y_fit) + c * np.cosh(d * y_fit) + e
        return p, y_fit, fy_fit, fy_guess
    return p

def fig_setup():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    return fig, ax

@jit(nopython=True)
def res_hyp(guesses, y, fy):
    a, b, c, d, e = guesses
    fy_calc = (a * np.sinh(b * y) + c * np.cosh(d * y)) * e * \
              np.sqrt((1 - y**2))
    err = np.sum(np.abs(fy - fy_calc))
    return err

@jit(nopython=True)
def eval_re(guesses, y):
    a, b, c, d, e = guesses
    fy_calc = (a * np.exp(b * y) + c * np.exp(d * y)) * e * \
              np.sqrt((1 - y**2))
    return fy_calc

@jit(nopython=True)
def res_exp(guesses, y, fy):
    fy_calc = eval_re(guesses, y)
    err = np.sum(np.abs(fy - fy_calc))
    return err

@jit(nopython=True)
def eval_rec(guesses, y):
    a, b, c, d, e, f = guesses
    fy_calc = (a * np.exp(b * y) + c * np.exp(d * y) + f) * e * \
              np.sqrt((1 - y**2))
    return fy_calc

@jit(nopython=True)
def res_exp_con(guesses, y, fy):
    fy_calc = eval_rec(guesses, y)
    err = np.sum(np.abs(fy - fy_calc))
    return err

@jit(nopython=True)
def eval_rel(guesses, y):
    a, b, c, d, e, f, g = guesses
    fy_calc = (a * np.exp(b * y) + c * np.exp(d * y) + f + g * y) * e * \
              np.sqrt((1 - y**2))
    return fy_calc

@jit(nopython=True)
def res_exp_lin(guesses, y, fy):
    fy_calc = eval_rel(guesses, y)
    err = np.sum(np.abs(fy - fy_calc))
    return err

@jit(nopython=True)
def eval_rhc(guesses, y):
    a, b, c, d, e, f = guesses
    fy_calc = (a * np.sinh(b * y) + c * np.cosh(d * y) + f) * e * \
              np.sqrt((1 - y**2))
    return fy_calc

@jit(nopython=True)
def res_hyp_con(guesses, y, fy):
    fy_calc = eval_rhc(guesses, y)
    err = np.sum(np.abs(fy - fy_calc))
    return err

@jit(nopython=True)
def eval_rhl(guesses, y):
    a, b, c, d, e, f, g = guesses
    fy_calc = (a * np.sinh(b * y) + c * np.cosh(d * y) + f + g * y) * e * \
              np.sqrt((1 - y**2))
    return fy_calc

@jit(nopython=True)
def res_hyp_lin(guesses, y, fy):
    fy_calc = eval_rhl(guesses, y)
    err = np.sum(np.abs(fy - fy_calc))
    return err

#try doing differential evolution fits:
def fit_data(bounds, y, fy, f_eval, f_err, y_fit=None):
    optres = opt.differential_evolution(f_err, bounds, args=(y, fy))
    if y_fit is not None:
        fy_fit = f_eval(optres.x, y_fit)
        return optres, fy_fit
    return optres

#Now try putting into polar coordinates ####doesn't help
def get_polar(y, fy):
    """Puts y, fy into polar coordinates"""
    dists = np.sqrt(y**2 + fy**2)
    angles = np.arctan2(fy, y)
    return dists, angles

y_fit = np.linspace(-1, 1, 1000)
max_b = 1e8
bounds_rhc = [(-max_b, max_b), (-10, 10), (-max_b, max_b), (-10, 10),
              (-max_b, max_b), (-max_b, max_b)]
bounds_rel = bounds_rhc + [(-max_b, max_b)]
th_i, mu_i = 0, 50

#fit_attempts = [fit_data(bounds_rhc[:-1], yvals, fy[th_i, mu_i, :], eval_re,
#                         res_exp, y_fit) for i in range(10)]
#fit_attempts = [fit_data(bounds_rhc, yvals, fy[th_i, mu_i, :], eval_rec,
#                         res_exp_con, y_fit) for i in range(10)]
#fit_attempts = [fit_data(bounds_rel, yvals, fy[th_i, mu_i, :], eval_rel,
#                         res_exp_lin, y_fit) for i in range(10)]
#fit_attempts = [fit_data(bounds_rel, yvals, fy[th_i, mu_i, :], eval_rhl,
#                         res_hyp_lin, y_fit) for i in range(10)]
#
#fig, ax = fig_setup()
#ax.plot(yvals, fy[th_i, mu_i, :], 'o')
#for i in range(len(fit_attempts)):
#    ax.plot(y_fit, fit_attempts[i][1])

#fit_attempts1 = [fit_data(bounds_rhc, yvals, fy[th_i, mu_i, :], eval_rhc,
#                         res_hyp_con, y_fit) for i in range(10)]
#fig1, ax1, = fig_setup()
#ax1.plot(yvals, fy[th_i, mu_i, :], 'o')
#for i in range(len(fit_attempts1)):
#    ax1.plot(y_fit, fit_attempts1[i][1])

#eval_rhc and res_hyp_con seem to work ok, so now to churn through a muval
