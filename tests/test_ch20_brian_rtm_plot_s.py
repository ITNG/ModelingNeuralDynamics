from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_rtm_plot_s_matches_matlab():
    """make_figure.m runs tau_r=0.2 then tau_r=1, reusing the same s/v
    arrays -- workspace at the end holds the tau_r=1 run.
    """
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    sm = ns.simulate_RTM_neuron_s(1.0 * ns.b2.uA, 1 * ns.b2.ms, 100 * ns.b2.ms)

    ref = run_matlab_script(ROOT / "matlab" / "20" / "RTM_PLOT_S", "make_figure.m", ["t", "v", "s"])

    rmse_v = trace_rmse(ref["t"], ref["v"], sm.t / ns.b2.ms, sm.vm[0] / ns.b2.mV)
    rmse_s = trace_rmse(ref["t"], ref["s"], sm.t / ns.b2.ms, sm.s[0])
    assert rmse_v < 2.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f} mV"
    assert rmse_s < 0.01, f"s RMSE vs MATLAB too high: {rmse_s:.4f}"
