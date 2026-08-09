from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script, spike_times

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_m_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter09.ipynb")
    import neurodynex3.tools.input_factory as input_factory

    current = input_factory.get_step_current(0, 300, ns.b2.ms, 1.5 * ns.b2.uA)
    sm = ns.simulate_RTM_M_neuron(current, 300 * ns.b2.ms)
    t = sm.t / ns.b2.ms
    v = sm.vm[0] / ns.b2.mV

    ref = run_matlab_script(ROOT / "matlab" / "09" / "RTM_M", "make_figure.m", ["t", "v"])

    spikes_test = spike_times(t, v)
    spikes_ref = spike_times(ref["t"], ref["v"])

    assert len(spikes_test) == len(spikes_ref) == 8
    np.testing.assert_allclose(spikes_test, spikes_ref, atol=2.0)
