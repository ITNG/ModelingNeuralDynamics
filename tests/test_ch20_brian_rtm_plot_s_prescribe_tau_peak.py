from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, load_python_port, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_brian_prescribed_tau_peak_matches_verified_python():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    py = load_python_port(
        ROOT / "python" / "20_Chemical_Synapses"
        / "RTM_PLOT_S_PRESCRIBE_TAU_PEAK" / "main.py"
    )
    tau_d_q = ns.release_time_constant(300.0, 10.0, 20.0)
    expected = py.tau_d_q_function(300.0, 10.0, 20.0)
    assert np.isclose(tau_d_q, expected, rtol=1e-8)
    assert np.isclose(ns.synaptic_peak_time(300.0, 10.0, tau_d_q), 20.0, atol=0.02)

    py.tau_d, py.tau_r, py.tau_d_q = 300.0, 10.0, expected
    py.i_ext = 0.12
    py.t = np.arange(0.0, 500.0, 0.01)
    py_trace = py.odeint(py.derivative, py.initial_condition(-70.0), py.t)
    sm = ns.simulate_RTM_neuron_q_s(
        0.12 * ns.b2.uA,
        tau_r=10 * ns.b2.ms,
        tau_d=300 * ns.b2.ms,
        tau_d_q=tau_d_q * ns.b2.ms,
        simulation_time=500 * ns.b2.ms,
    )
    assert trace_rmse(py.t, py_trace[:, 4], sm.t / ns.b2.ms, sm.s[0]) < 0.03
