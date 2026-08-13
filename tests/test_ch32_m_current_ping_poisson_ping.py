from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def _ns():
    return load_notebook_definitions_as_module(ROOT / "python" / "chapter32.ipynb")


def test_plot_phi_structure():
    ns = _ns()
    x, phi_vec = ns.simulate_plot_phi()
    assert np.all(phi_vec > 0) and np.all(phi_vec < 1)
    # phi(x) lies above the diagonal for small x and crosses it once
    assert phi_vec[0] > x[0]


def test_plot_psi_structure():
    ns = _ns()
    x, psi_vec = ns.simulate_plot_psi()
    assert np.all(psi_vec > 0) and np.all(psi_vec < 1)
    # psi is decreasing (an inhibitory PRC-like map)
    assert psi_vec[0] > psi_vec[-1]


def test_plot_psi_phi_structure():
    ns = _ns()
    x, psi_vec, phi_vec = ns.simulate_plot_psi_phi()
    assert np.all(psi_vec >= 0) and np.all(psi_vec <= 1)
    assert np.all(phi_vec >= 0) and np.all(phi_vec <= 1)


def test_m_current_ping_1_structure():
    ns = _ns()
    # the simulation is an explicit call, not an import side effect:
    # load_notebook_definitions_as_module only pulls imports/defs out of
    # the notebook's code cells, never top-level calls.
    assert not hasattr(ns, "t_e_spikes")
    assert not hasattr(ns, "v_plot")
    t_e, i_e, t_i, i_i, v_plot = ns.simulate_m_current_ping_1(t_final_run=200.0)
    # exact values below hold because simulate_m_current_ping_1 seeds its
    # own fresh rng internally, so every call (not just the first) is
    # independently reproducible.
    assert len(t_e) == 368
    assert len(t_i) == 282
    assert len(t_e) == len(i_e)
    assert len(t_i) == len(i_i)
    assert np.isclose(v_plot.max(), 48.377, rtol=0.0, atol=0.01)


def _check_ping(t_e, i_e, t_i, i_i, num_e=200, num_i=50, min_e_spikes=100, min_i_spikes=50):
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial E and I
    # firing) rather than exact spike times.
    assert len(t_e) > min_e_spikes
    assert len(t_i) > min_i_spikes


def test_m_current_ping_1_closeup_structure():
    ns = _ns()
    t_e, i_e, t_i, i_i, v_plot = ns.simulate_m_current_ping_closeup()
    _check_ping(t_e, i_e, t_i, i_i, min_e_spikes=500, min_i_spikes=1000)


def test_m_current_ping_2_closeup_structure():
    # the legacy _1/_2/_3_CLOSEUP scripts were parameter-identical, so the
    # notebook exposes a single simulate_m_current_ping_closeup() for all
    # three; call it again to represent the "second" closeup.
    ns = _ns()
    t_e, i_e, t_i, i_i, v_plot = ns.simulate_m_current_ping_closeup()
    _check_ping(t_e, i_e, t_i, i_i, min_e_spikes=500, min_i_spikes=1000)


def test_m_current_ping_3_closeup_structure():
    ns = _ns()
    t_e, i_e, t_i, i_i, v_plot = ns.simulate_m_current_ping_closeup()
    _check_ping(t_e, i_e, t_i, i_i, min_e_spikes=500, min_i_spikes=1000)


def test_m_current_ping_1_from_rest_structure():
    ns = _ns()
    t_e, i_e, t_i, i_i = ns.simulate_m_current_ping_1_from_rest()
    assert len(t_e) > 500
    assert len(t_i) > 1000
    # network is silent for the first ~200 ms (zero drive, from rest)
    assert t_e.min() > 1.0


def test_ping_clusters_structure():
    ns = _ns()
    t_e, i_e, t_i, i_i, v_plot = ns.simulate_ping_clusters()
    assert len(t_e) > 500
    assert len(t_i) > 500


def test_poisson_ping_1_structure():
    ns = _ns()
    t_e, i_e, t_i, i_i, v_one, num_e, num_i, m_steps = ns.simulate_poisson_ping_1()
    assert len(t_e) > 500
    assert len(t_i) > 100
    assert num_i == num_e // 4


def test_poisson_ping_2_structure():
    ns = _ns()
    t_e, i_e, t_i, i_i, v_one, num_e, num_i, m_steps = ns.simulate_poisson_ping_2()
    assert len(t_e) > 300
    assert len(t_i) > 300


def test_poisson_ping_3_structure():
    ns = _ns()
    t_e, i_e, t_i, i_i, v_one, num_e, num_i, m_steps = ns.simulate_poisson_ping_3()
    assert len(t_e) > 50
    assert len(t_i) > 300
    # strong, all-to-all I-I/E-I coupling with low heterogeneity produces
    # tightly synchronized I-cell volleys
    counts_i = np.bincount(i_i.astype(int), minlength=num_i)
    assert np.count_nonzero(counts_i) > 0.9 * num_i


def test_poisson_ping_3_voltage_trace_structure():
    ns = _ns()
    t_e, i_e, t_i, i_i, v_one, num_e, num_i, m_steps = ns.simulate_poisson_ping_3(track_cell=2)
    assert len(t_e) > 50
    assert len(t_i) > 300
    assert len(v_one) == m_steps + 1
    # the tracked cell spikes at least once (voltage crosses well above rest)
    assert v_one.max() > 0
