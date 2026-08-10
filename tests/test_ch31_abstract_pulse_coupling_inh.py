from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "31_ING_Rhythms/ABSTRACT_PULSE_COUPLING_INH"
MATLAB_DIR = "31/ABSTRACT_PULSE_COUPLING_INH"


def test_abstract_pulse_coupling_inh_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi", "ind"], timeout=30)

    assert np.allclose(py.phi, ref["phi"], atol=1e-9)
    # matlab's ind is 1-indexed into a 1000-length array; python's ind is
    # 0-indexed into the same array, so they should differ by exactly 1.
    assert np.array_equal(py.ind + 1, ref["ind"].astype(int))
