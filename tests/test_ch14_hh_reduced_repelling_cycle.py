from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "14_Model_Neurons_of_Bifurcation_Type_2/HH_REDUCED_REPELLING_CYCLE"
MATLAB_DIR = "14/HH_REDUCED_REPELLING_CYCLE"


def test_hh_reduced_repelling_cycle_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    # matlab's v/n at the end of the script are the reverse (repelling
    # cycle) trace, its last computed trajectory
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["v_c", "n_c", "v", "n"])

    assert np.isclose(py.v_c, ref["v_c"], rtol=1e-8)
    assert np.isclose(py.n_c, ref["n_c"], rtol=1e-8)

    t_ref = np.arange(len(ref["v"])) * py.dt
    t_py = np.arange(len(py.v_rep)) * py.dt
    rmse_v = trace_rmse(t_ref, ref["v"], t_py, py.v_rep)
    rmse_n = trace_rmse(t_ref, ref["n"], t_py, py.n_rep)
    assert rmse_v < 1.0, f"v RMSE vs MATLAB too high: {rmse_v:.4f}"
    assert rmse_n < 0.01, f"n RMSE vs MATLAB too high: {rmse_n:.5f}"
