from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script, spike_times
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_hh_refractoriness_matches_matlab():
    """make_figure.m's workspace at the end holds the last (onset=9) panel,
    same as tests/test_ch04_hh_refractoriness.py for the Python port.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter04.ipynb")
    t = ns.stm2.t / ns.b2.ms
    v = ns.stm2.vm[0] / ns.b2.mV

    ref = run_matlab_script(
        ROOT / "matlab" / "04" / "HH_REFRACTORINESS", "make_figure.m", ["t", "v"]
    )

    spikes_test = spike_times(t, v)
    spikes_ref = spike_times(ref["t"], ref["v"])

    assert len(spikes_test) == len(spikes_ref) == 4
    np.testing.assert_allclose(spikes_test, spikes_ref, atol=1.0)
