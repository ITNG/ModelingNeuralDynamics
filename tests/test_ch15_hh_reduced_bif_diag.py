from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "15_Canard_Explosions/HH_REDUCED_BIF_DIAG"
MATLAB_DIR = "15/HH_REDUCED_BIF_DIAG"


@pytest.mark.slow
def test_hh_reduced_bif_diag_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
        ["red_i", "red_v", "green_i", "green_v",
         "i_ext_uc", "Bmax_vec", "Bmin_vec",
         "i_ext_sc", "Amax_vec", "Amin_vec"],
        timeout=60,
    )

    red_i, red_v = zip(*py.red)
    assert np.allclose(red_i, ref["red_i"], atol=1e-6)
    assert np.allclose(red_v, ref["red_v"], atol=1e-6)

    green_i, green_v = zip(*py.green)
    assert np.allclose(green_i, ref["green_i"], atol=1e-6)
    assert np.allclose(green_v, ref["green_v"], atol=1e-6)

    i_sc, amax, amin = zip(*py.sc)
    assert np.allclose(i_sc, ref["i_ext_sc"], atol=1e-6)
    assert np.allclose(amax, ref["Amax_vec"], atol=1e-6)
    assert np.allclose(amin, ref["Amin_vec"], atol=1e-6)

    # unstable-cycle branch: only the (few) points where the backward
    # integration stays bounded on both sides -- the "if Bmax<200" filter
    # matlab and python both apply. matlab's saved i_ext_uc/Bmax_vec/
    # Bmin_vec already have the plotting script's prepended first point
    # (the stable-cycle amplitude at i_ext_sc[0]), which main.py only adds
    # in its __main__ block -- prepend the same way here for comparison.
    i_uc, bmax, bmin = zip(*py.uc)
    i_uc = (i_sc[0], *i_uc)
    bmax = (amax[0], *bmax)
    bmin = (amin[0], *bmin)
    assert np.allclose(i_uc, ref["i_ext_uc"], atol=1e-6)
    assert np.allclose(bmax, ref["Bmax_vec"], atol=1e-3)
    assert np.allclose(bmin, ref["Bmin_vec"], atol=1e-3)
