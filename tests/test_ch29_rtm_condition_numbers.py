from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "29_Stability_of_the_Synchronous_State/RTM_CONDITION_NUMBERS"
MATLAB_DIR = "29/RTM_CONDITION_NUMBERS"


def test_rtm_condition_numbers_matches_matlab():
    # matlab's script computes P_base/percentage_increase twice (weak/slow
    # then strong/fast synapse), reusing the variable names, so only the
    # second (strong/fast) block's values survive in the final workspace.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_table.m",
                             ["P_base", "P_lowered_I", "percentage_increase"], timeout=60)

    assert np.isclose(py.P_base_strong, ref["P_base"], atol=0.5)
    assert np.isclose(py.pct_strong['i_ext'], ref["percentage_increase"], atol=0.2)
