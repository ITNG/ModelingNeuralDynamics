from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "33_M_Current_PING_and_PINB"


def test_m_current_beta_with_gj_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing)
    # rather than exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "M_CURRENT_BETA_WITH_GJ" / "main.py")
    assert len(py.t_e_spikes) > 200
    assert py.num_e == 200


def _check_ping(name, min_e_spikes, min_i_spikes, num_e=40, num_i=10):
    py = load_python_port(ROOT / "python" / PYTHON_BASE / name / "main.py")
    assert len(py.t_e_spikes) > min_e_spikes
    assert len(py.t_i_spikes) > min_i_spikes
    assert py.num_e == num_e
    assert py.num_i == num_i


def test_m_current_ping_4_structure():
    _check_ping("M_CURRENT_PING_4", min_e_spikes=50, min_i_spikes=30)


def test_m_current_ping_5_structure():
    _check_ping("M_CURRENT_PING_5", min_e_spikes=50, min_i_spikes=30)


def test_m_current_ping_6_structure():
    _check_ping("M_CURRENT_PING_6", min_e_spikes=100, min_i_spikes=100)


def test_m_current_ping_7_structure():
    _check_ping("M_CURRENT_PING_7", min_e_spikes=100, min_i_spikes=100)


def test_m_current_ping_8_structure():
    _check_ping("M_CURRENT_PING_8", min_e_spikes=100, min_i_spikes=100)


def _check_pinb(name):
    py = load_python_port(ROOT / "python" / PYTHON_BASE / name / "main.py")
    assert len(py.t_e_spikes) > 500
    assert len(py.t_i_spikes) > 100
    assert py.num_e == 200
    assert py.num_i == 50
    assert len(py.lfp) == py.m_steps + 1


def test_pinb_1_structure():
    _check_pinb("PINB_1")


def test_pinb_2_structure():
    _check_pinb("PINB_2")


def test_pinb_3_structure():
    _check_pinb("PINB_3")
