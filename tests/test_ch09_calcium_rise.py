from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_calcium_rise_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    v, r = ns.calcium_rise()

    ref = run_matlab_script(ROOT / "matlab" / "09/CALCIUM_RISE", "make_figure.m", ["r"])

    assert np.allclose(r, ref["r"], rtol=1e-6, atol=1e-6)
