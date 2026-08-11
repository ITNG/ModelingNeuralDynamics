from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "14/HH_REDUCED_COUNT_FP"


@pytest.mark.slow
def test_hh_reduced_count_fp_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter14.ipynb")
    num_fp, i_ext_vec = ns.simulate_hh_reduced_count_fp()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "main.m", ["num_fp"])

    assert num_fp.max() == int(ref["num_fp"].max())
    assert num_fp.min() == int(ref["num_fp"].min())
