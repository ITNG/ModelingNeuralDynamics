from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "21/LIF_NETWORK_WITH_GJ"


@pytest.mark.slow
def test_lif_network_with_gj_matches_matlab():
    # matlab reuses v for both runs (coupled, then uncoupled) -- only the
    # second (uncoupled, epsilon=0) trace survives to the end of the script
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter21.ipynb")
    v_uncoupled, spikes_uncoupled = ns.simulate_lif_network_with_gj(0.0)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(v_uncoupled, ref["v"], atol=1e-6)
