from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "17_Frequency_Current_Curves/ERISIR_F_I_CURVE"
MATLAB_DIR = "17/ERISIR_F_I_CURVE"


@pytest.mark.slow
def test_erisir_f_i_curve_matches_matlab():
    # ~80s each side: forward+backward sweep over 31 i_ext values, each
    # settling (or spiking 4x) from the previous i_ext's final state.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["f_vec", "I_c", "I_star"], timeout=180)

    # matlab's f_vec is overwritten in place across all three sweeps
    # (forward, then backward over the same array) -- final state is the
    # backward sweep result
    assert np.allclose(py.f_backward, ref["f_vec"], atol=1e-2)
    assert np.isclose(py.I_c, ref["I_c"], atol=1e-6)
    assert np.isclose(py.I_star, ref["I_star"], atol=1e-6)
