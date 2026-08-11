from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "24/RTM_SYNC"


@pytest.mark.slow
def test_rtm_sync_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter24.ipynb")
    t_spikes, i_spikes = ns.simulate_rtm_sync()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t_spikes", "i_spikes"])

    assert np.allclose(t_spikes, ref["t_spikes"], atol=1e-6)
    assert np.array_equal(i_spikes, ref["i_spikes"].astype(int))
