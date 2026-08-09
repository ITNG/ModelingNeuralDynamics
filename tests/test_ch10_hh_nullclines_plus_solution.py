from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "10_The_Slow_Fast_Phase_Plane/HH_NULLCLINES_PLUS_SOLUTION"
MATLAB_DIR = "10/HH_NULLCLINES_PLUS_SOLUTION"


def test_hh_nullclines_plus_solution_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["n_vec", "v", "n"])

    # n_vec: bisection-solved nullcline, same algorithm both sides -> tight.
    assert np.allclose(py.n_vec, ref["n_vec"], atol=1e-6)

    # v: odeint vs matlab's Heun/RK2 -- loose, same reasoning as elsewhere.
    t_ref = np.arange(len(ref["v"])) * py.dt
    rmse_v = trace_rmse(t_ref, ref["v"], py.t, py.v)
    assert rmse_v < 5.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
