from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script, spike_times
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_hh_voltage_trace_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter01.ipynb")
    t = ns.state_monitor.t / ns.b2.ms
    v = ns.state_monitor.vm[0] / ns.b2.mV

    ref = run_matlab_script(
        ROOT / "matlab" / "01" / "HH_VOLTAGE_TRACE", "make_figure.m", ["t", "v"]
    )

    spikes_test = spike_times(t, v)
    spikes_ref = spike_times(ref["t"], ref["v"])

    assert len(spikes_test) == len(spikes_ref) == 12
    np.testing.assert_allclose(spikes_test, spikes_ref, atol=1.5)
