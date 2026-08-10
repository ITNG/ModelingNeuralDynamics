from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "16_Model_Neurons_of_Bifurcation_Type_3/INAPIK_FIXED_POINTS"
MATLAB_DIR = "16/INAPIK_FIXED_POINTS"


@pytest.mark.slow
def test_inapik_fixed_points_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["i_green", "v_green", "i_magenta", "v_magenta"])

    i_green, v_green = zip(*py.points['green'])
    assert np.allclose(i_green, ref["i_green"], atol=1e-6)
    assert np.allclose(v_green, ref["v_green"], atol=1e-6)

    i_magenta, v_magenta = zip(*py.points['magenta'])
    assert np.allclose(i_magenta, ref["i_magenta"], atol=1e-6)
    assert np.allclose(v_magenta, ref["v_magenta"], atol=1e-6)
