from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "31/ABSTRACT_PULSE_COUPLING_INH"


@pytest.mark.slow
def test_abstract_pulse_coupling_inh_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter31.ipynb")
    phi, ind = ns.simulate_abstract_pulse_coupling_inh()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["phi", "ind"], timeout=30)

    assert np.allclose(phi, ref["phi"], atol=1e-9)
    # matlab's ind is 1-indexed into a 1000-length array; python's ind is
    # 0-indexed into the same array, so they should differ by exactly 1.
    assert np.array_equal(ind + 1, ref["ind"].astype(int))
