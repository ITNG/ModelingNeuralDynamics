from pathlib import Path

import numpy as np

from matlab_ref import run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "13/HOPF_SUB"


def test_hopf_sub_matches_matlab():
    r = np.arange(101) / 100 * 1.2
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["r"])
    assert np.allclose(r, ref["r"], atol=1e-9)
