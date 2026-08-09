from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "15_Canard_Explosions/FITZHUGH_NAGUMO_MICRO"
MATLAB_DIR = "15/FITZHUGH_NAGUMO_MICRO"


def test_fitzhugh_nagumo_micro_matches_matlab():
    # ~8 minutes in python (a full 10s trajectory for every i_ext where the
    # fixed point is an unstable spiral, to trace the canard-cycle
    # envelope), a couple minutes in matlab -- expensive but this is the
    # whole point of the "canard explosion" zoom.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
        ["i_red", "v_c_red", "i_green", "v_c_green",
         "i_cycle", "max_cycle", "min_cycle"],
        timeout=300,
    )

    i_red, v_red = zip(*py.branches['red'])
    assert np.allclose(i_red, ref["i_red"], atol=1e-6)
    assert np.allclose(v_red, ref["v_c_red"], atol=1e-6)

    i_green, v_green = zip(*py.branches['green'])
    assert np.allclose(i_green, ref["i_green"], atol=1e-6)
    assert np.allclose(v_green, ref["v_c_green"], atol=1e-6)

    i_cycle, vmax, vmin = zip(*py.cycle)
    assert np.allclose(i_cycle, ref["i_cycle"], atol=1e-6)
    assert np.allclose(vmax, ref["max_cycle"], atol=1e-3)
    assert np.allclose(vmin, ref["min_cycle"], atol=1e-3)
