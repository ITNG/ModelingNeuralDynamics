from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "21_Gap_Junctions/WB_NETWORK_WITH_GJ_SUBTHRESHOLD"
MATLAB_DIR = "21/WB_NETWORK_WITH_GJ_SUBTHRESHOLD"


@pytest.mark.slow
def test_wb_network_with_gj_subthreshold_matches_matlab():
    # matlab reuses v for both runs (always-on gap junction, then
    # subthreshold-only) -- only the second trace survives to the end
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(py.v_sub, ref["v"], atol=1e-4)
