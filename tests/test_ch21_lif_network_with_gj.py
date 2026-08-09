from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "21_Gap_Junctions/LIF_NETWORK_WITH_GJ"
MATLAB_DIR = "21/LIF_NETWORK_WITH_GJ"


def test_lif_network_with_gj_matches_matlab():
    # matlab reuses v for both runs (coupled, then uncoupled) -- only the
    # second (uncoupled, epsilon=0) trace survives to the end of the script
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    assert np.allclose(py.v_uncoupled, ref["v"], atol=1e-6)
