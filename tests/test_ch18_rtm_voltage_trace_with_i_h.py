from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "18_Bistability_Resulting_from_Rebound_Firing/RTM_VOLTAGE_TRACE_WITH_I_H"
MATLAB_DIR = "18/RTM_VOLTAGE_TRACE_WITH_I_H"


def test_rtm_voltage_trace_with_i_h_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    t = np.arange(len(py.v)) * py.dt
    rmse = trace_rmse(t, ref["v"], t, py.v)
    assert rmse < 1.0, f"v RMSE vs MATLAB too high: {rmse:.4f}"
