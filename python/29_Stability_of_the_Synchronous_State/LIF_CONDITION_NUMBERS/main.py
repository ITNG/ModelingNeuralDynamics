import numpy as np

dt = 0.01
dt05 = dt / 2


def P0(tau_m, J, g, tau_I):
    I = J + 1 / tau_m
    v = 0.
    k = 1
    while v < 1:
        v_old = v
        v_inc = -v / tau_m + I - g * np.exp(-(k - 1) * dt / tau_I) * v
        v_tmp = v + dt05 * v_inc
        v_inc = -v_tmp / tau_m + I - g * np.exp(-(k - 0.5) * dt / tau_I) * v_tmp
        v = v + dt * v_inc
        k += 1
    period = ((k - 2) * dt * (v - 1) + (k - 1) * dt * (1 - v_old)) / (v - v_old)
    return period


def P1(tau_m, J, g, tau_I):
    I = J + 1 / tau_m
    if g <= J:
        return 0.
    v = 1.
    k = 1
    while v <= 1:
        v_old = v
        v_inc = -v / tau_m + I - g * np.exp(-(k - 1) * dt / tau_I) * v
        v_tmp = v + dt05 * v_inc
        v_inc = -v_tmp / tau_m + I - g * np.exp(-(k - 0.5) * dt / tau_I) * v_tmp
        v = v + dt * v_inc
        k += 1
    period = ((k - 2) * dt * (v - 1) + (k - 1) * dt * (1 - v_old)) / (v - v_old)
    return period


def P(tau_m, J, g, tau_I):
    return (P0(tau_m, J, g, tau_I) + P1(tau_m, J, g, tau_I)) / 2


tau_m = 10.

results = {}
for label, (J, g, tau_I) in (
    ('slow', (0.02, 0.15, 9.)),
    ('fast', (0.02, 2., 1.)),
):
    P_here = P(tau_m, J, g, tau_I)
    pct_J = (P(tau_m, 0.99 * J, g, tau_I) - P_here) / P_here * 100
    pct_g = (P(tau_m, J, g * 1.01, tau_I) - P_here) / P_here * 100
    pct_tau_I = (P(tau_m, J, g, tau_I * 1.01) - P_here) / P_here * 100
    results[label] = {'P': P_here, 'J': pct_J, 'g': pct_g, 'tau_I': pct_tau_I}

if __name__ == "__main__":
    for label, r in results.items():
        print(label, r)
