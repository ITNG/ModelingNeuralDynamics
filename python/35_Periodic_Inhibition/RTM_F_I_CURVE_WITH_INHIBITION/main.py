import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
from numba import njit
from numba.extending import register_jitable

c = 1.
g_k, g_na, g_l = 80., 100., 0.1
v_k, v_na, v_l = -100., 50., -67.

alpha = 1.
Period = 25.
g_bar = 0.1  # strength of tonic inhibition
N = 2000
denom = np.mean(np.exp(alpha * np.cos(np.pi * np.arange(N) / N) ** 2) - 1)
v_rev = -75.  # reversal potential of tonic/periodic inhibition

i_ext_low, i_ext_high = 0., 3.
i_ext_vec = i_ext_low + np.arange(101) / 100 * (i_ext_high - i_ext_low)

dt = 0.01
dt05 = dt / 2


@register_jitable
def m_inf(v):
    alpha_m = 0.32 * (v + 54) / (1 - exp(-(v + 54) / 4))
    beta_m = 0.28 * (v + 27) / (exp((v + 27) / 5) - 1)
    return alpha_m / (alpha_m + beta_m)


@register_jitable
def alpha_h(v):
    return 0.128 * exp(-(v + 50) / 18)


@register_jitable
def beta_h(v):
    return 4. / (1 + exp(-(v + 27) / 5))


@register_jitable
def alpha_n(v):
    return 0.032 * (v + 52) / (1 - exp(-(v + 52) / 5))


@register_jitable
def beta_n(v):
    return 0.5 * exp(-(v + 57) / 40)


@register_jitable
def g_periodic(t, amplitude=g_bar):
    return amplitude * (np.exp(alpha * np.cos(np.pi * t / Period) ** 2) - 1) / denom


@register_jitable
def step(v, m, h, n, i_ext, g_inhib_old, g_inhib_mid):
    v_inc = (g_na * m ** 3 * h * (v_na - v) + g_k * n ** 4 * (v_k - v) + g_l * (v_l - v)
              + g_inhib_old * (v_rev - v) + i_ext) / c
    h_inc = alpha_h(v) * (1 - h) - beta_h(v) * h
    n_inc = alpha_n(v) * (1 - n) - beta_n(v) * n

    v_tmp = v + dt05 * v_inc
    m_tmp = m_inf(v_tmp)
    h_tmp = h + dt05 * h_inc
    n_tmp = n + dt05 * n_inc

    v_inc = (g_na * m_tmp ** 3 * h_tmp * (v_na - v_tmp) + g_k * n_tmp ** 4 * (v_k - v_tmp)
              + g_l * (v_l - v_tmp) + g_inhib_mid * (v_rev - v_tmp) + i_ext) / c
    h_inc = alpha_h(v_tmp) * (1 - h_tmp) - beta_h(v_tmp) * h_tmp
    n_inc = alpha_n(v_tmp) * (1 - n_tmp) - beta_n(v_tmp) * n_tmp

    v_new = v + dt * v_inc
    m_new = m_inf(v_new)
    h_new = h + dt * h_inc
    n_new = n + dt * n_inc
    return v_new, m_new, h_new, n_new


def _f_i_curve_tonic_python(i_ext_values):
    '''For each drive in i_ext_values (continuing from the final state of
    the previous drive), integrate until either a fixed point is
    detected (frequency 0) or the 4th spike occurs (frequency from the
    3rd-to-4th interspike interval).'''
    N_check = round(1000 / dt)
    v, m, h, n = -70., m_inf(-70.), 0.7, 0.6
    f_vec = np.zeros(len(i_ext_values))

    for ijk, i_ext in enumerate(i_ext_values):
        v_hist = [v]
        m_hist = [m]
        h_hist = [h]
        n_hist = [n]
        t_spikes = []
        num_spikes = 0
        k = 0
        f = 0.
        while True:
            k += 1
            v, m, h, n = step(v_hist[-1], m_hist[-1], h_hist[-1], n_hist[-1], i_ext, g_bar, g_bar)
            v_hist.append(v)
            m_hist.append(m)
            h_hist.append(h)
            n_hist.append(n)

            if (k - 1) % N_check == 0 and k > 1:
                vv = np.array(v_hist[-N_check - 1:])
                mm = np.array(m_hist[-N_check - 1:])
                hh = np.array(h_hist[-N_check - 1:])
                nn = np.array(n_hist[-N_check - 1:])
                if (vv.max() - vv.min() < 1e-4 * abs(vv.max() + vv.min())
                        and mm.max() - mm.min() < 1e-4 * abs(mm.max() + mm.min())
                        and hh.max() - hh.min() < 1e-4 * abs(hh.max() + hh.min())
                        and nn.max() - nn.min() < 1e-4 * abs(nn.max() + nn.min())):
                    f = 0.
                    break

            if v_hist[-1] < -20 and v_hist[-2] >= -20:
                num_spikes += 1
                t_spikes.append((k * dt * (20 + v_hist[-2]) + (k - 1) * dt * (-20 - v_hist[-1]))
                                 / (v_hist[-2] - v_hist[-1]))
            if num_spikes == 4:
                f = 1000. / (t_spikes[3] - t_spikes[2])
                break

        f_vec[ijk] = f
        v, m, h, n = v_hist[-1], m_hist[-1], h_hist[-1], n_hist[-1]

    return f_vec


def _f_i_curve_periodic_python(i_ext_values, g_amplitude):
    t_final = 2000.
    m_steps = round(t_final / dt)
    t = np.arange(m_steps + 1) * dt
    g_store = g_periodic(t, g_amplitude)

    v, m, h, n = -70., m_inf(-70.), 0.7, 0.6
    f_vec = np.zeros(len(i_ext_values))

    for ijk, i_ext in enumerate(i_ext_values):
        num_spikes = 0
        v_old = v
        for k in range(1, m_steps + 1):
            v, m, h, n = step(v, m, h, n, i_ext, g_store[k - 1], (g_store[k - 1] + g_store[k]) / 2)
            if v < -20 and v_old >= -20:
                num_spikes += 1
            v_old = v
        f_vec[ijk] = num_spikes / t_final * 1000

    return f_vec


# njit without cache=True: these scripts are exec'd from a file path (not
# imported as real modules), so numba's on-disk cache can't re-import them
# to rebuild the compiled environment. Compilation is redone per process.
_f_i_curve_tonic_jit = njit(_f_i_curve_tonic_python)
_f_i_curve_periodic_jit = njit(_f_i_curve_periodic_python)


def compute_f_i_curves(i_ext_values=i_ext_vec, use_numba=True):
    '''F-I curves under tonic and under periodic inhibition, one entry
    per drive in i_ext_values.'''
    if use_numba:
        tonic_kernel = _f_i_curve_tonic_jit
        periodic_kernel = _f_i_curve_periodic_jit
    else:
        tonic_kernel = _f_i_curve_tonic_python
        periodic_kernel = _f_i_curve_periodic_python
    return tonic_kernel(i_ext_values), periodic_kernel(i_ext_values, g_bar)


if __name__ == "__main__":
    f_vec_tonic, f_vec_periodic = compute_f_i_curves()

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(i_ext_vec, f_vec_tonic, '-r', linewidth=2, label='tonic inhibition')
    ax.plot(i_ext_vec, f_vec_periodic, '-b', linewidth=2, label='periodic inhibition')
    ax.axis([i_ext_low, i_ext_high, 0, 100])
    ax.set_xlabel(r'$I$ [$\mu$A/cm$^2$]')
    ax.set_ylabel('$f$')
    ax.legend()

    plt.tight_layout()
    plt.savefig("fig.png")
    # plt.show()
