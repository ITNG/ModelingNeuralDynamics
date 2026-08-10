from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "04_Numerical_Solution_of_HH_ODEs/HH_LIMIT_CYCLE"
MATLAB_DIR = "04/HH_LIMIT_CYCLE"


@pytest.mark.slow
def test_hh_limit_cycle_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    t_p = py.np.arange(0, py.t_final, py.dt)
    sol = odeint(py.derivative, py.x0, t_p)
    v_p, n_p = sol[:, 0], sol[:, 2]

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v", "n"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    rmse_n = trace_rmse(ref["t"], ref["n"], t_p, n_p)
    assert rmse_v < 1.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f} mV"
    assert rmse_n < 0.01, f"n RMSE vs MATLAB too high: {rmse_n:.4f}"
