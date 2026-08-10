from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "12_Two_Dimensional_Bifurcation_Analysis/RTM_2D_INVARIANT_CYCLE"
MATLAB_DIR = "12/RTM_2D_INVARIANT_CYCLE"


@pytest.mark.slow
def test_rtm_2d_invariant_cycle_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    def f(v):
        return py.g_k * py.n_inf(v) ** 4 * (py.v_k - v) \
            + py.g_na * py.m_inf(v) ** 3 * py.h_inf(v) * (py.v_na - v) \
            + py.g_l * (py.v_l - v)

    v_star = py.find_root(f, -63, -62, 1e-14)
    v_0 = py.find_root(f, -75, -65, 1e-12)

    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v_star", "v_0", "v", "n"]
    )
    # closed-form-ish bisection roots, same algorithm both sides -> tight
    assert np.isclose(v_star, ref["v_star"], rtol=1e-9)
    assert np.isclose(v_0, ref["v_0"], rtol=1e-9)

    # matlab reuses the same v/n names throughout; the last assignment
    # (and so what's left in the workspace at the end) is the final sim
    # in panel 3 -- the one seeded at v_star-0.5 (the red trajectory)
    n_star = py.n_inf(v_star)
    v, n = py.simulate(v_star - 0.5, n_star - 0.005, i_ext=0., t_final=300.)
    t_ref = np.arange(len(ref["v"])) * 0.005
    t_py = np.arange(len(v)) * 0.005
    rmse = trace_rmse(t_ref, ref["v"], t_py, v)
    assert rmse < 5.0, f"v RMSE vs MATLAB too high: {rmse:.2f}"
