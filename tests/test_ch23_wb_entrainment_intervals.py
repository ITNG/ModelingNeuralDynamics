from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "23/WB_ENTRAINMENT_INTERVALS"


@pytest.mark.slow
def test_wb_entrainment_intervals_matches_matlab():
    # very slow on the MATLAB side (book's own comment: "This takes quite
    # a while to run") -- 201 g_syn values x 1e6-step WB-neuron sims. The
    # Python side is numba-jitted and fast; only the MATLAB reference call
    # needs the long timeout.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter23.ipynb")
    g_syn_vec, n_vec, sigma_vec = ns.simulate_wb_entrainment_intervals()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["g_syn_vec", "n_vec", "sigma_vec"], timeout=3600)

    assert np.allclose(g_syn_vec, ref["g_syn_vec"], atol=1e-9)
    assert np.allclose(n_vec, ref["n_vec"], atol=0.1)
    ind_py = sigma_vec < 1e-3
    ind_ref = ref["sigma_vec"] < 1e-3
    assert np.array_equal(ind_py, ind_ref)
