from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "16/INAPIK_FIXED_POINTS"


@pytest.mark.slow
def test_inapik_fixed_points_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter16.ipynb")
    points, i_ext_vec = ns.simulate_inapik_fixed_points()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["i_green", "v_green", "i_magenta", "v_magenta"])

    i_green, v_green = zip(*points['green'])
    assert np.allclose(i_green, ref["i_green"], atol=1e-6)
    assert np.allclose(v_green, ref["v_green"], atol=1e-6)

    i_magenta, v_magenta = zip(*points['magenta'])
    assert np.allclose(i_magenta, ref["i_magenta"], atol=1e-6)
    assert np.allclose(v_magenta, ref["v_magenta"], atol=1e-6)
