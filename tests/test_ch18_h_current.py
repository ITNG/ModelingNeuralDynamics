import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from matlab_ref import MATLAB, load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "18_Bistability_Resulting_from_Rebound_Firing/H_CURRENT"
MATLAB_DIR = ROOT / "matlab" / "18/H_CURRENT"


def test_h_current_matches_matlab():
    import os
    import shutil
    if not (os.path.exists(MATLAB) or shutil.which("matlab")):
        pytest.skip("MATLAB not available on this machine")

    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    # r_inf.m/tau_r.m are only ever called inline inside a plot() call in
    # make_figure.m, so evaluate them directly instead of running the script
    with tempfile.TemporaryDirectory() as tmp:
        outfile = Path(tmp) / "out.mat"
        cmd = (
            f"cd('{MATLAB_DIR}'); v=[-100:50]; "
            f"r=r_inf(v); t=tau_r(v); save('{outfile}', 'v', 'r', 't');"
        )
        subprocess.run([MATLAB, "-batch", cmd], check=True, timeout=60, capture_output=True, text=True)
        data = loadmat(outfile)

    v = data["v"].squeeze()
    assert np.allclose(py.v, v)
    assert np.allclose(py.r_inf(v), data["r"].squeeze(), atol=1e-9)
    assert np.allclose(py.tau_r(v), data["t"].squeeze(), atol=1e-9)
