from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "29_Stability_of_the_Synchronous_State/LIF_CONDITION_NUMBERS"
MATLAB_DIR = "29/LIF_CONDITION_NUMBERS"


@pytest.mark.slow
def test_lif_condition_numbers_matches_matlab():
    # matlab's script only prints unnamed "ans" values (semicolon-free
    # expressions), so only the very last one ("ans") survives in the
    # workspace to be saved -- the tau_I percentage of the second
    # (fast) block.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_table.m",
                             ["ans"], timeout=30)

    assert np.isclose(py.results['fast']['tau_I'], ref["ans"], atol=1e-3)
