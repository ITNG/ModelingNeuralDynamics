from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "28/WEAKLY_COUPLED_1"


@pytest.mark.slow
def test_weakly_coupled_1_matches_matlab():
    # matlab's script reuses t_vec/psi_vec/psi for both epsilon runs, so
    # by the end of the script they hold only the second (epsilon=0.1)
    # run's results.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter28.ipynb")
    t_vec_2, psi_vec_2, t_de_2, psi_de_2 = ns.simulate_weakly_coupled_1(epsilon=0.1)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t_vec", "psi_vec", "psi"], timeout=30)

    assert np.allclose(t_vec_2, ref["t_vec"], atol=1e-6)
    assert np.allclose(psi_vec_2, ref["psi_vec"], atol=1e-6)
    assert np.allclose(psi_de_2, ref["psi"], atol=1e-6)
