import numpy as np
import pylab as pl
from numpy import exp
from main import *


def alpha_h_o(v):
    return 0.07 * exp(-(v + 63.0) / 20.0)


def alpha_m_o(v):
    q = (v+38.0)/10.0
    return q / (1.0 - exp(-q))


def alpha_n_o(v):
    return 0.018 * (v - 25.0) / (1.0 - exp(-(v - 25.0) / 25.0))


def beta_h_o(v):
    return 1.0 / (exp(-(v + 33.0) / 10.0) + 1.0)


def beta_m_o(v):
    return 4.0 * exp(-(v + 65.0) / 18.0)


def beta_n_o(v):
    return 0.0036 * (35.0 - v) / (1.0 - exp(-(35.0 - v) / 12.0))


def b_o_inf(v):
    return 1.0/(1.0+exp((v+71.0)/7.3))


def a_o_inf(v):
    return 1.0 / (1.0 + exp(-(v + 14.0) / 16.6))


def h_o_inf(v):
    return alpha_h_o(v) / (alpha_h_o(v) + beta_h_o(v))


def m_o_inf(v):
    return alpha_m_o(v) / (alpha_m_o(v) + beta_m_o(v))


def n_o_inf(v):
    return alpha_n_o(v) / (alpha_n_o(v) + beta_n_o(v))


def r_o_inf(v):
    return 1.0 / (1.0 + exp((v + 84.0) / 10.2))


def tau_a_o(v):
    return 5.0


def tau_b_o(v):
    return 1. / (0.000009 / exp((v - 26) / 28.5) +
                 0.014 / (0.2 + exp(-(v + 70.0) / 11.0)))


def tau_r_o(v):
    return 1.0 / (exp(-14.59 - 0.086 * v) + exp(-1.87 + 0.0701 * v))


def derivative(x0, t):

    v, h, n, r, a, b = x0

    I_K = g_k * n ** 4 * (v - v_k)
    I_Na = g_na * h * m_o_inf(v) ** 3 * (v - v_na)
    I_L = g_l * (v - v_l)
    I_H = g_h * r * (v - v_h)
    I_A = g_A * a * b * (v - v_A)

    dv = (i_ext - I_L - I_K - I_Na - I_H - I_A) / c
    dh = alpha_h_o(v) * (1.0 - h) - beta_h_o(v) * h
    dn = alpha_n_o(v) * (1.0 - n) - beta_n_o(v) * n
    dr = (r_o_inf(v) - r) / tau_r_o(v)
    da = (a_o_inf(v) - a) / tau_a_o(v)
    db = (b_o_inf(v) - b) / tau_b_o(v)

    return [dv, dh, dn, dr, da, db]
