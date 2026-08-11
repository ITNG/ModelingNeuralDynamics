from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "15/FITZHUGH_NAGUMO_MACRO"


@pytest.mark.slow
def test_fitzhugh_nagumo_macro_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter15.ipynb")
    branches, i_c, cycle_pts, i_ext_vec = ns.simulate_fitzhugh_nagumo_macro()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["i_red", "v_c_red", "i_green", "v_c_green", "i_c"],
                             timeout=60)

    # for these parameters (a=5, tau_n=60) the fixed point is always a
    # spiral (stable=red or unstable=green), never a node/saddle -- black/
    # blue/magenta branches are empty on both sides
    assert branches['black'] == [] and branches['blue'] == [] \
        and branches['magenta'] == []

    i_red, v_red = zip(*branches['red'])
    assert np.allclose(i_red, ref["i_red"], atol=1e-6)
    assert np.allclose(v_red, ref["v_c_red"], atol=1e-6)

    i_green, v_green = zip(*branches['green'])
    assert np.allclose(i_green, ref["i_green"], atol=1e-6)
    assert np.allclose(v_green, ref["v_c_green"], atol=1e-6)

    # the I where the fixed point loses stability (last crossing found in
    # the loop -- matlab overwrites i_c each time it fires)
    assert np.isclose(i_c, ref["i_c"], atol=1e-6)
