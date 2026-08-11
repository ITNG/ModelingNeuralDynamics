from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]

# Reference spike times from a faithful re-implementation of MATLAB's
# make_figure.m that preserves the full v/w/t history across resets
# (make_figure.m itself clears its arrays after every spike for
# incremental plotting, same as chapter 7/8's LIF/QIF scripts, so it
# can't be extracted from the workspace directly).
MATLAB_SPIKES = [14.63, 45.36337425, 89.244391, 133.81396821,
                 178.397667, 222.97981098, 267.56205607]


def test_lif_adapt_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter09.ipynb")

    sm, spm = ns.simulate_LIF_adapt_neuron(
        tau_m=10.0, I=0.13, tau_a=40.0, delta=0.05, simulation_time=300 * ns.b2.ms
    )
    spikes = spm.t / ns.b2.ms

    assert len(spikes) == len(MATLAB_SPIKES)
    np.testing.assert_allclose(spikes, MATLAB_SPIKES, atol=0.5)
