from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "22/WILSON_COWAN_LOWERING_W_EE"


@pytest.mark.slow
def test_wilson_cowan_lowering_w_ee_matches_matlab():
    # matlab reuses E/I for all three w_EE runs -- only the last (w_EE=1.0)
    # trace survives to the end of the script
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter22.ipynb")
    E, I = ns.simulate_wilson_cowan_lowering_w_ee(1.0)
    m_steps = len(E) - 1
    ind = slice(round(4 * m_steps / 5) - 1, m_steps + 1)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["E", "I"])

    assert np.allclose(E[ind], ref["E"], atol=1e-6)
    assert np.allclose(I[ind], ref["I"], atol=1e-6)
