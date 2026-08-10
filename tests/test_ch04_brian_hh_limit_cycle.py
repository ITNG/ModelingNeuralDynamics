from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_hh_limit_cycle_matches_matlab():
    """Same simulation as HH_SOLUTION, plotted as a v-n phase portrait.

    n doesn't have the huge dv/dt spike upstrokes that make trace_rmse
    unreliable for v, so it's a clean whole-trace check here.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter04.ipynb")
    t = ns.state_monitor.t / ns.b2.ms
    n = ns.state_monitor.n[0]

    ref = run_matlab_script(ROOT / "matlab" / "04" / "HH_LIMIT_CYCLE", "make_figure.m", ["t", "n"])

    rmse = trace_rmse(ref["t"], ref["n"], t, n)
    assert rmse < 0.05, f"n RMSE vs MATLAB too high: {rmse:.4f}"
