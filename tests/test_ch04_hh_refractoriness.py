from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "04_Numerical_Solution_of_HH_ODEs/HH_REFRACTORINESS"
MATLAB_DIR = "04/HH_REFRACTORINESS"


def test_hh_refractoriness_matches_matlab():
    """make_figure.m runs all three pulse-onset panels (-2, 5, 9) in one
    script, reusing v/t each time -- so its workspace at the end only holds
    the last (onset=9) trace. Compare that one against the matching Python
    run.
    """
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    t_p = py.np.arange(0, py.t_final, py.dt)
    v_p = odeint(py.derivative, py.x0, t_p, args=(10.0, 9.0))[:, 0]

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"])

    rmse = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    assert rmse < 1.0, f"RMSE vs MATLAB too high: {rmse:.2f} mV"
