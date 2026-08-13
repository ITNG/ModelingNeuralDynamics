from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "python" / "chapter33.ipynb"


def _ns():
    return load_notebook_definitions_as_module(NOTEBOOK)


def test_m_current_beta_with_gj_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing)
    # rather than exact spike times.
    ns = _ns()
    t_e_spikes, i_e_spikes, num_e = ns.simulate_m_current_beta_with_gj()
    assert len(t_e_spikes) > 200
    assert num_e == 200


def _check_ping(fn, min_e_spikes, min_i_spikes, num_e=40, num_i=10):
    t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes = fn()
    assert len(t_e_spikes) > min_e_spikes
    assert len(t_i_spikes) > min_i_spikes
    assert i_e_spikes.max() < num_e
    assert i_i_spikes.max() < num_i


def test_m_current_ping_4_structure():
    ns = _ns()
    _check_ping(ns.simulate_m_current_ping_4, min_e_spikes=50, min_i_spikes=30)


def test_m_current_ping_5_structure():
    ns = _ns()
    _check_ping(ns.simulate_m_current_ping_5, min_e_spikes=50, min_i_spikes=30)


def test_m_current_ping_6_structure():
    ns = _ns()
    _check_ping(ns.simulate_m_current_ping_6, min_e_spikes=100, min_i_spikes=100)


def test_m_current_ping_7_structure():
    ns = _ns()
    _check_ping(ns.simulate_m_current_ping_7, min_e_spikes=100, min_i_spikes=100)


def test_m_current_ping_8_structure():
    ns = _ns()
    _check_ping(ns.simulate_m_current_ping_8, min_e_spikes=100, min_i_spikes=100)


def _check_pinb(fn):
    t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, lfp, m_steps = fn()
    assert len(t_e_spikes) > 500
    assert len(t_i_spikes) > 100
    assert i_e_spikes.max() < 200
    assert i_i_spikes.max() < 50
    assert len(lfp) == m_steps + 1


def test_pinb_1_structure():
    ns = _ns()
    _check_pinb(ns.simulate_pinb_1)


def test_pinb_2_structure():
    ns = _ns()
    _check_pinb(ns.simulate_pinb_2)


def test_pinb_3_structure():
    ns = _ns()
    _check_pinb(ns.simulate_pinb_3)
