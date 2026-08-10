from pathlib import Path

import numpy as np

from matlab_ref import (
    load_notebook_definitions_as_module,
    load_python_port,
    trace_rmse,
)

ROOT = Path(__file__).resolve().parents[1]


def test_brian_rtm_plot_q_matches_verified_python():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    py = load_python_port(
        ROOT / "python" / "20_Chemical_Synapses" / "RTM_PLOT_Q" / "main.py"
    )
    py.tau_r, py.tau_d = 0.1, 2.0
    py.t = np.arange(0.0, 100.0, 0.01)
    py_trace = py.odeint(py.derivative, py.initial_condition(-70.0), py.t)
    sm = ns.simulate_RTM_neuron_q_s(
        1.0 * ns.b2.uA,
        tau_r=0.1 * ns.b2.ms,
        tau_d=2.0 * ns.b2.ms,
        tau_d_q=2.0 * ns.b2.ms,
        simulation_time=100 * ns.b2.ms,
    )
    assert trace_rmse(py.t, py_trace[:, 0], sm.t / ns.b2.ms, sm.vm[0] / ns.b2.mV) < 2.0
    assert trace_rmse(py.t, py_trace[:, 3], sm.t / ns.b2.ms, sm.q[0]) < 0.03
    assert trace_rmse(py.t, py_trace[:, 4], sm.t / ns.b2.ms, sm.s[0]) < 0.03
