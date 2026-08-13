from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "29/ILLUSTRATE_P0_AND_P1"


@pytest.mark.slow
def test_illustrate_p0_and_p1_matches_matlab():
    # matlab's script reuses the variable name "period" for both runs, so
    # by the end of the script it holds only the second (P1) run's value.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter29.ipynb")
    _, _, _, _, _, period_1 = ns.simulate_p0_and_p1()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["period"], timeout=30)

    assert np.isclose(period_1, ref["period"], atol=1e-6)
