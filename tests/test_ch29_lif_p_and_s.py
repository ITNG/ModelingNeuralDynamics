from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "29/LIF_P_AND_S"


@pytest.mark.slow
def test_lif_p_and_s_matches_matlab():
    # matlab's script reuses P_vec/S_vec for all three sweeps, so by the
    # end of the script they hold only the third (varying J) sweep.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter29.ipynb")
    (_, _, _, _, _, _, J_vec, P_vec_J, S_vec_J) = ns.simulate_lif_p_and_s()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["J_vec", "P_vec", "S_vec"], timeout=60)

    assert np.allclose(J_vec, ref["J_vec"], atol=1e-9)
    assert np.allclose(P_vec_J, ref["P_vec"], atol=1e-2)
    assert np.allclose(S_vec_J, ref["S_vec"], atol=1e-3)
