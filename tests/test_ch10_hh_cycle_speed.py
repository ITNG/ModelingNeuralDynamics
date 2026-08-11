from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "10/HH_CYCLE_SPEED"


@pytest.mark.slow
def test_hh_cycle_speed_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter10.ipynb")
    v_vec, n_vec, v, n, speed, n_inf = ns.simulate_hh_cycle_speed()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["n_vec", "v", "n"])

    # n_vec: bisection-solved nullcline, same algorithm both sides -> tight.
    assert np.allclose(n_vec, ref["n_vec"], atol=1e-6)

    # v/n: Heun/RK2 trajectory, tail end of a 150ms sim at dt=0.001 -- loose,
    # same reasoning as other spiking-trace comparisons in this repo.
    dt = 0.001
    t_ref = np.arange(len(ref["v"])) * dt
    t_py = np.arange(len(v)) * dt
    rmse_v = trace_rmse(t_ref, ref["v"], t_py, v)
    assert rmse_v < 5.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
