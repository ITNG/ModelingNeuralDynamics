from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def _ns():
    return load_notebook_definitions_as_module(ROOT / "python" / "chapter37.ipynb")


def test_no_reset_structure():
    # deterministic (no RNG) -- exact match confirmed visually against
    # matlab's figure.pdf.
    ns = _ns()
    tau_m, tau_m_hat = ns.no_reset_time_constants()
    assert tau_m_hat < tau_m  # inhibition speeds up the membrane


def test_ping_thr_1_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing,
    # a monotonic drive ramp) rather than exact spike times.
    ns = _ns()
    t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, lfp, i_ext_e, num_e, num_i = ns.simulate_ping_thr_1()
    assert len(t_e_spikes) > 200
    assert len(t_i_spikes) > 5
    assert np.all(np.diff(i_ext_e) > 0)  # deterministic ramp


def test_ping_thr_1_zoom_structure():
    ns = _ns()
    t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, lfp, i_ext_e, num_e, num_i = ns.simulate_ping_thr_1_zoom(
        t_final=500.0)
    mask = (i_e_spikes >= 71) & (i_e_spikes <= 77)
    assert mask.sum() > 20


def test_thresholding_structure():
    # deterministic bisection search (no RNG). Reference values below are
    # matlab's own hardcoded output, copied verbatim from ratios.m in the
    # matlab source (the first 5 of its 6 entries correspond to the 5
    # g_bar values swept by make_table.m); our independently computed w
    # matches to 5-6 significant digits.
    ns = _ns()
    g_bar_vec, w, t = ns.threshold_width_sweep()
    w_matlab = np.array([0.036266326904297, 0.012943267822266, 0.004600524902344,
                          0.001636505126953, 5.893707275390625e-04])
    assert np.allclose(w, w_matlab, rtol=1e-3)
    # w decreases roughly geometrically with g_bar (ratio ~2.8 each step)
    ratios = w[:4] / w[1:]
    assert np.all((ratios > 2.5) & (ratios < 3.1))
