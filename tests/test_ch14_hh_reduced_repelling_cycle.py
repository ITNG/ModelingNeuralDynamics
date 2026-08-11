from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "14/HH_REDUCED_REPELLING_CYCLE"


@pytest.mark.slow
def test_hh_reduced_repelling_cycle_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter14.ipynb")
    v_c, n_c, v_attr, n_attr, v_rep, n_rep, dt = ns.simulate_hh_reduced_repelling_cycle()

    # matlab's v/n at the end of the script are the reverse (repelling
    # cycle) trace, its last computed trajectory
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["v_c", "n_c", "v", "n"])

    assert np.isclose(v_c, ref["v_c"], rtol=1e-8)
    assert np.isclose(n_c, ref["n_c"], rtol=1e-8)

    t_ref = np.arange(len(ref["v"])) * dt
    t_py = np.arange(len(v_rep)) * dt
    rmse_v = trace_rmse(t_ref, ref["v"], t_py, v_rep)
    rmse_n = trace_rmse(t_ref, ref["n"], t_py, n_rep)
    assert rmse_v < 1.0, f"v RMSE vs MATLAB too high: {rmse_v:.4f}"
    assert rmse_n < 0.01, f"n RMSE vs MATLAB too high: {rmse_n:.5f}"
