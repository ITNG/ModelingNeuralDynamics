import numpy as np
import pylab as pl
from numpy import exp
from main import *

def alpha_h(v):
    return 0.07 * exp(-(v + 63.0) / 20.0)


def alpha_m(v):
    q = (v+38.0)/10.0
    return q / (1.0 - exp(-q))


def alpha_n(v):
    return 0.018 * (v - 25.0) / (1.0 - exp(-(v - 25.0) / 25.0))


def beta_h(v):
    return 1.0 / (exp(-(v + 33.0) / 10.0) + 1.0)


def beta_m(v):
    return 4.0 * exp(-(v + 65.0) / 18.0)


def beta_n(v):
    return 0.0036 * (35.0 - v) / (1.0 - exp(-(35.0 - v) / 12.0))


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))


def tau_h(v):
    return 1.0 / (alpha_h(v) + beta_h(v))


def tau_m(v):
    return 1.0 / (alpha_m(v) + beta_m(v))


def tau_n(v):
    return 1.0 / (alpha_n(v) + beta_n(v))


def derivative(x0, t):
    v, h, n = x0
    
    dv = i_ext - g_na * h * m_inf(v) ** 3 * \
        (v - v_na) - g_k * n ** 4 * (v - v_k) - g_l * (v - v_l)
    dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n
    dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h

    return [dv, dh, dn]
