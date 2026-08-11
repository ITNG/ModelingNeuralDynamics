from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "21/WB_NETWORK_WITH_GJ_SUBTHRESHOLD"


@pytest.mark.slow
def test_wb_network_with_gj_subthreshold_matches_matlab():
    # matlab reuses v for both runs (always-on gap junction, then
    # subthreshold-only) -- only the second trace survives to the end
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter21.ipynb")
    gap_always_on = np.array([[0.0, 0.01], [0.01, 0.0]])

    def gap_subthreshold(v1, thr=-50.0):
        return gap_always_on if v1 <= thr else np.zeros((2, 2))

    v_sub = ns.simulate_wb_network_with_gj_subthreshold(gap_subthreshold)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(v_sub, ref["v"], atol=1e-4)
