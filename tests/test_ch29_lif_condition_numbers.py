from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "29/LIF_CONDITION_NUMBERS"


@pytest.mark.slow
def test_lif_condition_numbers_matches_matlab():
    # matlab's script only prints unnamed "ans" values (semicolon-free
    # expressions), so only the very last one ("ans") survives in the
    # workspace to be saved -- the tau_I percentage of the second
    # (fast) block.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter29.ipynb")
    results = ns.compute_lif_condition_numbers()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_table.m",
                             ["ans"], timeout=30)

    assert np.isclose(results['fast']['tau_I'], ref["ans"], atol=1e-3)
