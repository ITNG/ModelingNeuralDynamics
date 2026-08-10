from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_wb_network_with_gj_matches_matlab():
    py = load_python_port(ROOT / "python" / "21_Gap_Junctions" / "WB_NETWORK_WITH_GJ" / "main.py")
    t = py.np.arange(0, py.t_final, py.dt)
    sol = odeint(py.derivative, py.x0, t)
    v1, v2 = sol[:, 0], sol[:, 1]

    ref = run_matlab_script(ROOT / "matlab" / "21" / "WB_NETWORK_WITH_GJ", "make_figure.m", ["t", "v"])

    rmse_v1 = trace_rmse(ref["t"], ref["v"][0], t, v1)
    rmse_v2 = trace_rmse(ref["t"], ref["v"][1], t, v2)
    assert rmse_v1 < 1.5, f"v1 RMSE vs MATLAB too high: {rmse_v1:.2f} mV"
    assert rmse_v2 < 0.05, f"v2 RMSE vs MATLAB too high: {rmse_v2:.4f} mV"
