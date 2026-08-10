from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "28_Weakly_Coupled_Oscillators/WEAKLY_COUPLED_2"
MATLAB_DIR = "28/WEAKLY_COUPLED_2"


@pytest.mark.slow
def test_weakly_coupled_2_matches_matlab():
    # matlab's script reuses t_vec/psi_vec/psi for both epsilon runs, so
    # by the end of the script they hold only the second (epsilon=0.1)
    # run's results.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t_vec", "psi_vec", "psi"], timeout=30)

    assert np.allclose(py.t_vec_2, ref["t_vec"], atol=1e-6)
    assert np.allclose(py.psi_vec_2, ref["psi_vec"], atol=1e-6)
    assert np.allclose(py.psi_de_2, ref["psi"], atol=1e-6)
