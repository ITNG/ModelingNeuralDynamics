from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "28/WEAKLY_COUPLED_HETEROGENEOUS_1"


@pytest.mark.slow
def test_weakly_coupled_heterogeneous_1_matches_matlab():
    # matlab's script does "clear" between the c=0.08 and c=0.12 runs, so
    # only the second (c=0.12) run's t_vec/psi_vec/psi survive to the
    # final workspace.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter28.ipynb")
    t_de_2, psi_de_2, t_vec_2, psi_vec_2 = ns.simulate_weakly_coupled_heterogeneous_1(c=0.12)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t_vec", "psi_vec", "psi"], timeout=30)

    assert np.allclose(t_vec_2, ref["t_vec"], atol=1e-6)
    assert np.allclose(psi_vec_2, ref["psi_vec"], atol=1e-6)
    assert np.allclose(psi_de_2, ref["psi"], atol=1e-6)
