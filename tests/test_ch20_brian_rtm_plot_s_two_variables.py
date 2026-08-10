from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_plot_s_two_variables_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    sm = ns.simulate_RTM_neuron_q_s(0.12 * ns.b2.uA, tau_r=100 * ns.b2.ms, tau_d=300 * ns.b2.ms,
                                     tau_d_q=100 * ns.b2.ms, simulation_time=2000 * ns.b2.ms)

    ref = run_matlab_script(
        ROOT / "matlab" / "20" / "RTM_PLOT_S_TWO_VARIABLES", "make_figure.m", ["t", "v", "s"]
    )

    rmse_v = trace_rmse(ref["t"], ref["v"], sm.t / ns.b2.ms, sm.vm[0] / ns.b2.mV)
    rmse_s = trace_rmse(ref["t"], ref["s"], sm.t / ns.b2.ms, sm.s[0])
    assert rmse_v < 1.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f} mV"
    assert rmse_s < 0.01, f"s RMSE vs MATLAB too high: {rmse_s:.4f}"
