from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "14/HH_REDUCED_FP_EVS"


@pytest.mark.slow
def test_hh_reduced_fp_evs_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter14.ipynb")
    real_part, imag_part, i_ext_vec = ns.simulate_hh_reduced_fp_evs()
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["real_part", "imag_part"]
    )

    assert np.allclose(real_part, ref["real_part"], atol=1e-6)
    assert np.allclose(imag_part, ref["imag_part"], atol=1e-6)
