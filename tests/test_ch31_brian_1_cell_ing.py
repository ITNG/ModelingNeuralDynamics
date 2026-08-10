from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port, spike_times, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_brian_single_cell_ing_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    py = load_python_port(ROOT / "python" / "31_ING_Rhythms" / "1_CELL_ING" / "main.py")
    py.tau_dq_i = py.tau_d_q_function(py.tau_d_i, py.tau_r_i, py.tau_peak_i)
    t = np.arange(0.0, py.t_final, py.dt)
    sol = py.odeint(py.derivative, py.x0, t)
    expected_spikes = spike_times(t, sol[:, 0], threshold=-20.0)

    result = ns.simulate_single_cell_ing()
    actual_spikes = np.asarray(result["spike_ms"])

    assert trace_rmse(t, sol[:, 0], result["t_ms"], result["v_mv"]) < 3.0
    assert len(actual_spikes) == len(expected_spikes)
    assert np.isclose(result["period_ms"], np.diff(expected_spikes)[-1], atol=0.5)
