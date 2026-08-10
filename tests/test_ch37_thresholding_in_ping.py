from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "37_Thresholding_in_PING"


def test_no_reset_structure():
    # deterministic script (no RNG) -- exact match confirmed visually
    # against matlab's figure.pdf.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "NO_RESET" / "main.py")
    assert py.tau_m_hat < py.tau_m  # inhibition speeds up the membrane


def test_ping_thr_1_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing,
    # a monotonic drive ramp) rather than exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_THR_1" / "main.py")
    assert len(py.t_e_spikes) > 200
    assert len(py.t_i_spikes) > 5
    assert np.all(np.diff(py.i_ext_e) > 0)  # deterministic ramp


def test_ping_thr_1_zoom_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_THR_1_ZOOM" / "main.py")
    assert py.t_final == 500.0
    mask = (py.i_e_spikes >= 71) & (py.i_e_spikes <= 77)
    assert mask.sum() > 20


def test_thresholding_structure():
    # deterministic bisection search (no RNG). Reference values below are
    # matlab's own hardcoded output, copied verbatim from ratios.m in the
    # matlab source (the first 5 of its 6 entries correspond to the 5
    # g_bar values swept by make_table.m); our independently computed w
    # matches to 5-6 significant digits.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "THRESHOLDING" / "main.py")
    w_matlab = np.array([0.036266326904297, 0.012943267822266, 0.004600524902344,
                          0.001636505126953, 5.893707275390625e-04])
    assert np.allclose(py.w, w_matlab, rtol=1e-3)
    # w decreases roughly geometrically with g_bar (ratio ~2.8 each step)
    ratios = py.w[:4] / py.w[1:]
    assert np.all((ratios > 2.5) & (ratios < 3.1))
