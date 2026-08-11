from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "12/RTM_2D_INVARIANT_CYCLE"


@pytest.mark.slow
def test_rtm_2d_invariant_cycle_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter12.ipynb")

    v_star, v_0 = ns.rtm_2d_invariant_cycle_roots()

    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v_star", "v_0", "v", "n"]
    )
    # closed-form-ish bisection roots, same algorithm both sides -> tight
    assert np.isclose(v_star, ref["v_star"], rtol=1e-9)
    assert np.isclose(v_0, ref["v_0"], rtol=1e-9)

    # matlab reuses the same v/n names throughout; the last assignment
    # (and so what's left in the workspace at the end) is the final sim
    # in panel 3 -- the one seeded at v_star-0.5 (the red trajectory)
    n_star = ns.n_inf(v_star)
    v, n = ns.simulate_rtm_trajectory(v_star - 0.5, n_star - 0.005, i_ext=0.0, t_final=300.0)
    t_ref = np.arange(len(ref["v"])) * 0.005
    t_py = np.arange(len(v)) * 0.005
    rmse = trace_rmse(t_ref, ref["v"], t_py, v)
    assert rmse < 5.0, f"v RMSE vs MATLAB too high: {rmse:.2f}"
