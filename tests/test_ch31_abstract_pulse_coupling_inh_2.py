from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "31_ING_Rhythms/ABSTRACT_PULSE_COUPLING_INH_2"
MATLAB_DIR = "31/ABSTRACT_PULSE_COUPLING_INH_2"


def test_abstract_pulse_coupling_inh_2_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["varphi"], timeout=30)

    assert np.allclose(py.varphi, ref["varphi"], atol=1e-9)
    assert np.allclose(py.g(py.varphi), -py.epsilon * py.varphi * np.tanh((1 - py.varphi) * py.a) / np.tanh(py.a),
                        atol=1e-9)
