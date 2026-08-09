from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "05_The_Simple_Model_of_Neurons_in_Rodent_Brains/RTM_VOLTAGE_TRACE"
MATLAB_DIR = "05/RTM_VOLTAGE_TRACE"


def test_rtm_voltage_trace_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    t_p = py.np.arange(0, py.t_final, py.dt)
    v_p = odeint(py.derivative, py.x0, t_p)[:, 0]

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"])

    rmse = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    assert rmse < 2.0, f"RMSE vs MATLAB too high: {rmse:.2f} mV"
