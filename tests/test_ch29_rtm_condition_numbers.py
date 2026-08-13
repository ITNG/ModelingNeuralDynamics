from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "29/RTM_CONDITION_NUMBERS"


@pytest.mark.slow
def test_rtm_condition_numbers_matches_matlab():
    # matlab's script computes P_base/percentage_increase twice (weak/slow
    # then strong/fast synapse) AND three times within each block (i_ext,
    # g_syn, tau_syn perturbations), all reusing the same variable names --
    # so only the very last-computed value (strong/fast block, tau_syn
    # perturbation) survives in matlab's final workspace.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter29.ipynb")
    results = ns.compute_rtm_condition_numbers()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_table.m",
                             ["P_base", "P_lowered_I", "percentage_increase"], timeout=60)

    assert np.isclose(results['strong']['P_base'], ref["P_base"], atol=0.5)
    assert np.isclose(results['strong']['pct']['tau_syn'], ref["percentage_increase"], atol=0.2)
