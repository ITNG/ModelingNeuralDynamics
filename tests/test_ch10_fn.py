from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "10/FN"


@pytest.mark.slow
def test_fn_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter10.ipynb")
    dt = 0.01
    t_p, v_p, n_p = ns.simulate_fn(dt=dt)

    # matlab/10/FN/make_figure.m never names its time vector -- it's built
    # inline as `[0:m_steps]*dt` right in the plot call -- so only `v`
    # (the trajectory, overwriting the earlier nullcline sweep variable of
    # the same name) is available to pull from the workspace.
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])
    t_ref = np.arange(len(ref["v"])) * dt

    rmse = trace_rmse(t_ref, ref["v"], t_p, v_p)
    assert rmse < 1.0, f"RMSE vs MATLAB too high: {rmse:.2f}"
