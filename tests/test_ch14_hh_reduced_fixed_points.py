from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "14_Model_Neurons_of_Bifurcation_Type_2/HH_REDUCED_FIXED_POINTS"
MATLAB_DIR = "14/HH_REDUCED_FIXED_POINTS"


def test_hh_reduced_fixed_points_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
        ["v_c_vec", "real_part", "imag_part"],
    )

    v_c_vec = np.array([py.find_fixed_point(i) for i in py.i_ext_vec])
    assert np.allclose(v_c_vec, ref["v_c_vec"], rtol=1e-8)

    real_part = np.array([np.linalg.eigvals(py.jacobian(v, py.n_inf(v)))[0].real
                           for v in v_c_vec])
    assert np.allclose(real_part, ref["real_part"], atol=1e-6)
