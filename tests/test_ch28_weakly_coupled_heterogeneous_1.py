from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "28_Weakly_Coupled_Oscillators/WEAKLY_COUPLED_HETEROGENEOUS_1"
MATLAB_DIR = "28/WEAKLY_COUPLED_HETEROGENEOUS_1"


def test_weakly_coupled_heterogeneous_1_matches_matlab():
    # matlab's script does "clear" between the c=0.08 and c=0.12 runs, so
    # only the second (c=0.12) run's t_vec/psi_vec/psi survive to the
    # final workspace.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t_vec", "psi_vec", "psi"], timeout=30)

    assert np.allclose(py.t_vec_2, ref["t_vec"], atol=1e-6)
    assert np.allclose(py.psi_vec_2, ref["psi_vec"], atol=1e-6)
    assert np.allclose(py.psi_de_2, ref["psi"], atol=1e-6)
