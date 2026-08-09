from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "15_Canard_Explosions/FITZHUGH_NAGUMO_MACRO"
MATLAB_DIR = "15/FITZHUGH_NAGUMO_MACRO"


def test_fitzhugh_nagumo_macro_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["i_red", "v_c_red", "i_green", "v_c_green", "i_c"],
                             timeout=60)

    # for these parameters (a=5, tau_n=60) the fixed point is always a
    # spiral (stable=red or unstable=green), never a node/saddle -- black/
    # blue/magenta branches are empty on both sides
    assert py.branches['black'] == [] and py.branches['blue'] == [] \
        and py.branches['magenta'] == []

    i_red, v_red = zip(*py.branches['red'])
    assert np.allclose(i_red, ref["i_red"], atol=1e-6)
    assert np.allclose(v_red, ref["v_c_red"], atol=1e-6)

    i_green, v_green = zip(*py.branches['green'])
    assert np.allclose(i_green, ref["i_green"], atol=1e-6)
    assert np.allclose(v_green, ref["v_c_green"], atol=1e-6)

    # the I where the fixed point loses stability (last crossing found in
    # the loop -- matlab overwrites i_c each time it fires)
    assert np.isclose(py.i_c, ref["i_c"], atol=1e-6)
