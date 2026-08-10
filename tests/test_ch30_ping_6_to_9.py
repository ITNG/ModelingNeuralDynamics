from pathlib import Path

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "30_The_PING_Model_of_Gamma_Rhythms"


def test_ping_6_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing in
    # all three connectivity panels) rather than exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_6" / "main.py")
    for suffix in ("_1", "_2", "_3"):
        t_e = getattr(py, f"t_e_spikes{suffix}")
        t_i = getattr(py, f"t_i_spikes{suffix}")
        assert len(t_e) > 1000
        assert len(t_i) > 200


def test_ping_7_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_7" / "main.py")
    assert len(py.t_e_spikes) > 200
    assert len(py.t_i_spikes) > 30
    assert py.lfp.min() > -100 and py.lfp.max() < 50


def test_ping_8_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_8" / "main.py")
    assert len(py.t_e_spikes) > 200
    assert len(py.t_i_spikes) > 30
    assert py.lfp.min() > -100 and py.lfp.max() < 50


def test_ping_9_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_9" / "main.py")
    assert len(py.t_e_spikes) > 200
    assert len(py.t_i_spikes) > 30
    assert py.lfp.min() > -100 and py.lfp.max() < 50
