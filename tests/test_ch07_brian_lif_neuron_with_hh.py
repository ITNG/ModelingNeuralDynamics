from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script, spike_times
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_lif_neuron_with_hh_matches_matlab():
    """Same HH subthreshold trace underlies Figures 7.1-7.3 (LIF_NEURON_WITH_HH,
    TAU_M_FOR_HH, SUBTHR_FOR_HH); this checks the shared simulation's spike
    times, same repetitive-firing caveat as the ch01/ch04/ch05 brian tests.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter07.ipynb")

    ref = run_matlab_script(
        ROOT / "matlab" / "07" / "LIF_NEURON_WITH_HH", "make_figure.m", ["t", "v"]
    )

    spikes_test = spike_times(ns.t, ns.v, threshold=-20)
    spikes_ref = spike_times(ref["t"], ref["v"], threshold=-20)

    assert len(spikes_test) == len(spikes_ref) == 5
    np.testing.assert_allclose(spikes_test, spikes_ref, atol=1.0)
