from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "14/HH_REDUCED_FIXED_POINTS"


@pytest.mark.slow
def test_hh_reduced_fixed_points_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter14.ipynb")
    i_ext_vec = np.arange(1001) / 1000 * 15
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
        ["v_c_vec", "real_part", "imag_part"],
    )

    v_c_vec = np.array([ns.hh_reduced_find_fixed_point(i) for i in i_ext_vec])
    assert np.allclose(v_c_vec, ref["v_c_vec"], rtol=1e-8)

    real_part = np.array([np.linalg.eigvals(ns.hh_reduced_jacobian(v, ns.hh_reduced_n_inf(v)))[0].real
                           for v in v_c_vec])
    assert np.allclose(real_part, ref["real_part"], atol=1e-6)
