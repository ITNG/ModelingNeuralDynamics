from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script, spike_times, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_olm_with_h_current_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter34.ipynb")

    ref = run_matlab_script(
        ROOT / "matlab" / "34" / "OLM_WITH_H_CURRENT", "make_figure.m", ["t", "v", "r"]
    )

    rmse_r = trace_rmse(ref["t"], ref["r"], ns.sm_h.t / ns.b2.ms, ns.sm_h.r[0])
    assert rmse_r < 0.001, f"r RMSE vs MATLAB too high: {rmse_r:.5f}"

    spikes_test = spike_times(ns.sm_h.t / ns.b2.ms, ns.sm_h.vm[0] / ns.b2.mV, threshold=-20)
    spikes_ref = spike_times(ref["t"], ref["v"], threshold=-20)
    assert len(spikes_test) == len(spikes_ref) == 10
    np.testing.assert_allclose(spikes_test, spikes_ref, atol=0.5)
