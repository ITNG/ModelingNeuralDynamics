from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "10/HH_NULLCLINES_PLUS_SOLUTION"


@pytest.mark.slow
def test_hh_nullclines_plus_solution_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter10.ipynb")
    dt = 0.01
    v_vec, n_vec, t, v, n, sol, n_inf = ns.simulate_hh_nullclines_plus_solution(dt=dt)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["n_vec", "v", "n"])

    # n_vec: bisection-solved nullcline, same algorithm both sides -> tight.
    assert np.allclose(n_vec, ref["n_vec"], atol=1e-6)

    # v: odeint vs matlab's Heun/RK2 -- loose, same reasoning as elsewhere.
    t_ref = np.arange(len(ref["v"])) * dt
    rmse_v = trace_rmse(t_ref, ref["v"], t, v)
    assert rmse_v < 5.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
