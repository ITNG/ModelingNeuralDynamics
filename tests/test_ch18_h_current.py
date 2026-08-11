import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from matlab_ref import MATLAB, load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = ROOT / "matlab" / "18/H_CURRENT"


def test_h_current_matches_matlab():
    import os
    import shutil
    if not (os.path.exists(MATLAB) or shutil.which("matlab")):
        pytest.skip("MATLAB not available on this machine")

    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v, r, t = ns.simulate_h_current()

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

    v_ref = data["v"].squeeze()
    assert np.allclose(v, v_ref)
    assert np.allclose(ns.r_inf(v_ref), data["r"].squeeze(), atol=1e-9)
    assert np.allclose(ns.tau_r(v_ref), data["t"].squeeze(), atol=1e-9)
