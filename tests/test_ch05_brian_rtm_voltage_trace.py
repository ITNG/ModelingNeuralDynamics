from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script, spike_times
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_rtm_voltage_trace_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter05.ipynb")
    from mnd.brian.core import get_step_current

    current = get_step_current(0, 100, ns.b2.ms, 1.5 * ns.b2.uA)
    sm = ns.simulate_RTM_neuron(current, 100 * ns.b2.ms)
    t = sm.t / ns.b2.ms
    v = sm.vm[0] / ns.b2.mV

    ref = run_matlab_script(
        ROOT / "matlab" / "05" / "RTM_VOLTAGE_TRACE", "make_figure.m", ["t", "v"]
    )

    spikes_test = spike_times(t, v)
    spikes_ref = spike_times(ref["t"], ref["v"])

    assert len(spikes_test) == len(spikes_ref) == 6
    np.testing.assert_allclose(spikes_test, spikes_ref, atol=0.5)
