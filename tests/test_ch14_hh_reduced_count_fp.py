from pathlib import Path

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "14_Model_Neurons_of_Bifurcation_Type_2/HH_REDUCED_COUNT_FP"
MATLAB_DIR = "14/HH_REDUCED_COUNT_FP"


def test_hh_reduced_count_fp_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "main.m", ["num_fp"])

    assert py.num_fp.max() == int(ref["num_fp"].max())
    assert py.num_fp.min() == int(ref["num_fp"].min())
