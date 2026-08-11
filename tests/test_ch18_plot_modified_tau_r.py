import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from matlab_ref import MATLAB, load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = ROOT / "matlab" / "18/PLOT_MODIFIED_TAU_R"


def test_plot_modified_tau_r_matches_matlab():
    if not (os.path.exists(MATLAB) or shutil.which("matlab")):
        pytest.skip("MATLAB not available on this machine")

    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v, tau_original, tau_modified = ns.simulate_modified_tau_r()

    with tempfile.TemporaryDirectory() as tmp:
        outfile = Path(tmp) / "out.mat"
        cmd = (
            f"cd('{MATLAB_DIR}'); v=[-100:50]; "
            f"a=tau_r_original(v); b=tau_r_modified(v); save('{outfile}', 'v', 'a', 'b');"
        )
        subprocess.run([MATLAB, "-batch", cmd], check=True, timeout=60, capture_output=True, text=True)
        data = loadmat(outfile)

    v_ref = data["v"].squeeze()
    assert np.allclose(v, v_ref)
    assert np.allclose(tau_original, data["a"].squeeze(), atol=1e-9)
    assert np.allclose(tau_modified, data["b"].squeeze(), atol=1e-9)
