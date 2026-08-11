from pathlib import Path

from matlab_ref import (
    load_notebook_definitions_as_module,
    trace_rmse,
)

ROOT = Path(__file__).resolve().parents[1]


def test_brian_rtm_plot_q_matches_verified_python():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    ns_python = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    t, V, Q, S = ns_python.simulate_rtm_plot_q(tau_r=0.1, tau_d=2.0, t_final=100.0, dt=0.01)
    sm = ns.simulate_RTM_neuron_q_s(
        1.0 * ns.b2.uA,
        tau_r=0.1 * ns.b2.ms,
        tau_d=2.0 * ns.b2.ms,
        tau_d_q=2.0 * ns.b2.ms,
        simulation_time=100 * ns.b2.ms,
    )
    assert trace_rmse(t, V, sm.t / ns.b2.ms, sm.vm[0] / ns.b2.mV) < 2.0
    assert trace_rmse(t, Q, sm.t / ns.b2.ms, sm.q[0]) < 0.03
    assert trace_rmse(t, S, sm.t / ns.b2.ms, sm.s[0]) < 0.03
