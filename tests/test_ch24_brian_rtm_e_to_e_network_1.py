from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_e_to_e_network_1_matches_python_early_spikes():
    """python/24.../RTM_E_TO_E_NETWORK_1 starts 30 RTM cells in a splay
    state and shows them converge to a common (synchronized) firing
    pattern. Brian's rk4 (fixed dt) vs. that port's scipy.odeint (adaptive)
    integration diverge slowly cycle-to-cycle in this coupled 30-neuron
    network, so only the spike count and the first couple of spikes (before
    drift accumulates) are checked tightly.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter24.ipynb")
    b2 = ns.b2

    py_spikes_neuron0 = [0.83165858, 23.47236654, 48.01616887, 77.25025012,
                          109.05727813, 142.46542821, 176.94539525]

    brian_spikes_neuron0 = ns.spike_times_from_trace(ns.sm3.t / b2.ms, ns.sm3.vm[0] / b2.mV)

    assert len(brian_spikes_neuron0) == len(py_spikes_neuron0) == 7
    np.testing.assert_allclose(brian_spikes_neuron0[:2], py_spikes_neuron0[:2], atol=2.0)
