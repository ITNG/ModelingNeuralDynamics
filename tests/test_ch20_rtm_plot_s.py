from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_rtm_plot_s_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    t, V, S = ns.simulate_rtm_plot_s(tau_r=1.0, tau_d=2.0)

    ref = run_matlab_script(ROOT / "matlab" / "20" / "RTM_PLOT_S", "make_figure.m", ["t", "v", "s"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t, V)
    rmse_s = trace_rmse(ref["t"], ref["s"], t, S)
    assert rmse_v < 2.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f} mV"
    assert rmse_s < 0.01, f"s RMSE vs MATLAB too high: {rmse_s:.4f}"
