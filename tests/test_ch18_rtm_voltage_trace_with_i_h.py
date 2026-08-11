from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "18/RTM_VOLTAGE_TRACE_WITH_I_H"


@pytest.mark.slow
def test_rtm_voltage_trace_with_i_h_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v, t_final, dt = ns.simulate_rtm_voltage_trace_with_i_h()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    t = np.arange(len(v)) * dt
    rmse = trace_rmse(t, ref["v"], t, v)
    assert rmse < 1.0, f"v RMSE vs MATLAB too high: {rmse:.4f}"
