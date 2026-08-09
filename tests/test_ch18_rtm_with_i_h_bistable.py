from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "18_Bistability_Resulting_from_Rebound_Firing/RTM_WITH_I_H_BISTABLE"
MATLAB_DIR = "18/RTM_WITH_I_H_BISTABLE"


def test_rtm_with_i_h_bistable_matches_matlab():
    # matlab overwrites v/m/h/n/r for both runs -- only the second (firing)
    # trace survives to the end of the script
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "v_star"], timeout=60)

    assert np.isclose(py.v_star, ref["v_star"], atol=1e-6)
    t = np.arange(len(py.v_fire)) * py.dt
    rmse = trace_rmse(t, ref["v"], t, py.v_fire)
    assert rmse < 2.0, f"v RMSE vs MATLAB too high: {rmse:.4f}"
