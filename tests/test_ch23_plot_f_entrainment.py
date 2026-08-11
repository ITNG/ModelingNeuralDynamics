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
MATLAB_DIR = ROOT / "matlab" / "23/PLOT_F_ENTRAINMENT"


def test_plot_f_entrainment_matches_matlab():
    if not (os.path.exists(MATLAB) or shutil.which("matlab")):
        pytest.skip("MATLAB not available on this machine")

    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter23.ipynb")

    # F is a matlab function handle, only ever evaluated inline in
    # make_figure.m -- evaluate it directly instead of running the script
    with tempfile.TemporaryDirectory() as tmp:
        outfile = Path(tmp) / "out.mat"
        cmd = (
            f"cd('{MATLAB_DIR}'); T=25; tau=20; epsilon=0.2; "
            f"F=@(alpha) (alpha+epsilon)*exp(-T/tau); "
            f"a=F(0); b=F(1-epsilon); save('{outfile}', 'a', 'b');"
        )
        subprocess.run([MATLAB, "-batch", cmd], check=True, timeout=60, capture_output=True, text=True)
        data = loadmat(outfile)

    epsilon = 0.2
    assert np.isclose(ns.F_entrainment(0.0), data["a"].squeeze(), atol=1e-9)
    assert np.isclose(ns.F_entrainment(1 - epsilon), data["b"].squeeze(), atol=1e-9)
