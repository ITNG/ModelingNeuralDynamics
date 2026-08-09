from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "10_The_Slow_Fast_Phase_Plane/HH_CYCLE_SPEED"
MATLAB_DIR = "10/HH_CYCLE_SPEED"


def test_hh_cycle_speed_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["n_vec", "v", "n"])

    # n_vec: bisection-solved nullcline, same algorithm both sides -> tight.
    assert np.allclose(py.n_vec, ref["n_vec"], atol=1e-6)

    # v/n: Heun/RK2 trajectory, tail end of a 150ms sim at dt=0.001 -- loose,
    # same reasoning as other spiking-trace comparisons in this repo.
    t_ref = np.arange(len(ref["v"])) * py.dt
    t_py = np.arange(len(py.v)) * py.dt
    rmse_v = trace_rmse(t_ref, ref["v"], t_py, py.v)
    assert rmse_v < 5.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
