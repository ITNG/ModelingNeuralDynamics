import sys
from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "34_Nested_Gamma_Theta_Rhythms/PRE_OLM_VOLTAGE_TRACE"
MATLAB_DIR = "34/PRE_OLM_VOLTAGE_TRACE"


@pytest.mark.slow
def test_pre_olm_voltage_trace_matches_matlab():
    sys.path.insert(0, str(ROOT / "python" / PYTHON_DIR))
    try:
        py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    finally:
        sys.path.remove(str(ROOT / "python" / PYTHON_DIR))
    v0 = -63.0
    x0 = [v0, py.h_inf(v0), py.n_inf(v0)]
    t = py.np.arange(0, py.t_final, py.dt)
    v = odeint(py.derivative, x0, t)[:, 0]

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"])

    rmse = trace_rmse(ref["t"], ref["v"], t, v)
    assert rmse < 1.0, f"RMSE vs MATLAB too high: {rmse:.2f} mV"
