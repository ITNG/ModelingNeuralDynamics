from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "15/FITZHUGH_NAGUMO_MICRO"


@pytest.mark.slow
def test_fitzhugh_nagumo_micro_matches_matlab():
    # ~8 minutes in python (a full 10s trajectory for every i_ext where the
    # fixed point is an unstable spiral, to trace the canard-cycle
    # envelope), a couple minutes in matlab -- expensive but this is the
    # whole point of the "canard explosion" zoom.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter15.ipynb")
    branches, cycle, i_c_found, i_ext_vec = ns.simulate_fitzhugh_nagumo_micro()
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
        ["i_red", "v_c_red", "i_green", "v_c_green",
         "i_cycle", "max_cycle", "min_cycle"],
        timeout=300,
    )

    i_red, v_red = zip(*branches['red'])
    assert np.allclose(i_red, ref["i_red"], atol=1e-6)
    assert np.allclose(v_red, ref["v_c_red"], atol=1e-6)

    i_green, v_green = zip(*branches['green'])
    assert np.allclose(i_green, ref["i_green"], atol=1e-6)
    assert np.allclose(v_green, ref["v_c_green"], atol=1e-6)

    i_cycle, vmax, vmin = zip(*cycle)
    assert np.allclose(i_cycle, ref["i_cycle"], atol=1e-6)
    assert np.allclose(vmax, ref["max_cycle"], atol=1e-3)
    assert np.allclose(vmin, ref["min_cycle"], atol=1e-3)
