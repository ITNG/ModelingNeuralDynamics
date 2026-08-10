from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "23_Entrainment_by_Excitatory_Input_Pulses/WB_ENTRAINMENT_INTERVALS"
MATLAB_DIR = "23/WB_ENTRAINMENT_INTERVALS"


@pytest.mark.slow
def test_wb_entrainment_intervals_matches_matlab():
    # very slow both sides (book's own comment: "This takes quite a
    # while to run") -- 201 g_syn values x 1e6-step WB-neuron sims each.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["g_syn_vec", "n_vec", "sigma_vec"], timeout=3600)

    assert np.allclose(py.g_syn_vec, ref["g_syn_vec"], atol=1e-9)
    assert np.allclose(py.n_vec, ref["n_vec"], atol=0.1)
    ind_py = py.sigma_vec < 1e-3
    ind_ref = ref["sigma_vec"] < 1e-3
    assert np.array_equal(ind_py, ind_ref)
