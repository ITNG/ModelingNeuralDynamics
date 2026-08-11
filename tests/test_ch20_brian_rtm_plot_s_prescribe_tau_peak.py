from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_brian_prescribed_tau_peak_matches_verified_python():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    ns_python = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")

    tau_d_q = ns.release_time_constant(300.0, 10.0, 20.0)
    expected = ns_python.tau_d_q_function(tau_d=300.0, tau_r=10.0, tau_hat=20.0)
    assert np.isclose(tau_d_q, expected, rtol=1e-8)
    assert np.isclose(ns.synaptic_peak_time(300.0, 10.0, tau_d_q), 20.0, atol=0.02)

    t, V, Q, S = ns_python.simulate_two_stage_synapse(
        tau_r=10.0, tau_d=300.0, tau_d_q=expected, i_ext=0.12, t_final=500.0, dt=0.01)
    sm = ns.simulate_RTM_neuron_q_s(
        0.12 * ns.b2.uA,
        tau_r=10 * ns.b2.ms,
        tau_d=300 * ns.b2.ms,
        tau_d_q=tau_d_q * ns.b2.ms,
        simulation_time=500 * ns.b2.ms,
    )
    assert trace_rmse(t, S, sm.t / ns.b2.ms, sm.s[0]) < 0.03
