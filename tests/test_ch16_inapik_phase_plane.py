from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "16_Model_Neurons_of_Bifurcation_Type_3/INAPIK_PHASE_PLANE"
MATLAB_DIR = "16/INAPIK_PHASE_PLANE"


def test_inapik_phase_plane_matches_matlab():
    # matlab overwrites v/n/t_final/dt on every panel -- only the last
    # (I=4.4) trajectory survives to the end of the script
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "n"])

    p = py.panels[-1]
    t_ref = np.arange(len(ref["v"])) * p['dt']
    t_py = np.arange(len(p['v_full'])) * p['dt']
    rmse_v = trace_rmse(t_ref, ref["v"], t_py, p['v_full'])
    rmse_n = trace_rmse(t_ref, ref["n"], t_py, p['n_full'])
    assert rmse_v < 5.0, f"v RMSE vs MATLAB too high: {rmse_v:.3f}"
    assert rmse_n < 0.05, f"n RMSE vs MATLAB too high: {rmse_n:.4f}"
