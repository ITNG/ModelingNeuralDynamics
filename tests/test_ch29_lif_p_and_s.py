from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "29_Stability_of_the_Synchronous_State/LIF_P_AND_S"
MATLAB_DIR = "29/LIF_P_AND_S"


@pytest.mark.slow
def test_lif_p_and_s_matches_matlab():
    # matlab's script reuses P_vec/S_vec for all three sweeps, so by the
    # end of the script they hold only the third (varying J) sweep.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["J_vec", "P_vec", "S_vec"], timeout=60)

    assert np.allclose(py.J_vec, ref["J_vec"], atol=1e-9)
    assert np.allclose(py.P_vec_J, ref["P_vec"], atol=1e-2)
    assert np.allclose(py.S_vec_J, ref["S_vec"], atol=1e-3)
