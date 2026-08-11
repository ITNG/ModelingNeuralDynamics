from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_tau_d_q_function_reproduces_peak_time():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    tau_d_q = ns.tau_d_q_function(tau_d=300.0, tau_r=10.0, tau_hat=20.0)
    peak = ns.tau_peak_function(tau_d=300.0, tau_r=10.0, tau_d_q=tau_d_q)

    assert np.isclose(peak, 20.0, atol=1e-6)


def test_simulate_two_stage_synapse_with_prescribed_tau_d_q_smoke():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    tau_d_q = ns.tau_d_q_function(tau_d=300.0, tau_r=10.0, tau_hat=20.0)
    t, V, Q, S = ns.simulate_two_stage_synapse(tau_r=10.0, tau_d=300.0, tau_d_q=tau_d_q, t_final=50.0)

    assert np.all(np.isfinite(V)) and np.all(np.isfinite(S))
