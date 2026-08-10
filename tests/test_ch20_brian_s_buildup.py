from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_s_buildup_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter20.ipynb")
    sm = ns.simulate_RTM_neuron_q_s(0.2 * ns.b2.uA, tau_r=10 * ns.b2.ms, tau_d=300 * ns.b2.ms,
                                     tau_d_q=5 * ns.b2.ms, simulation_time=2000 * ns.b2.ms)

    ref = run_matlab_script(ROOT / "matlab" / "20" / "S_BUILDUP", "make_figure.m", ["t", "v", "s"])

    rmse_v = trace_rmse(ref["t"], ref["v"], sm.t / ns.b2.ms, sm.vm[0] / ns.b2.mV)
    rmse_s = trace_rmse(ref["t"], ref["s"], sm.t / ns.b2.ms, sm.s[0])
    assert rmse_v < 3.5, f"v RMSE vs MATLAB too high: {rmse_v:.2f} mV"
    assert rmse_s < 0.01, f"s RMSE vs MATLAB too high: {rmse_s:.4f}"
