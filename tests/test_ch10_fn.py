from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "10_The_Slow_Fast_Phase_Plane/FN"
MATLAB_DIR = "10/FN"


def test_fn_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    t_p = py.np.arange(0, py.t_final, py.dt)
    sol = odeint(py.derivative, py.x0, t_p)
    v_p = sol[:, 0]

    # matlab/10/FN/make_figure.m never names its time vector -- it's built
    # inline as `[0:m_steps]*dt` right in the plot call -- so only `v`
    # (the trajectory, overwriting the earlier nullcline sweep variable of
    # the same name) is available to pull from the workspace.
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])
    t_ref = py.np.arange(len(ref["v"])) * py.dt

    rmse = trace_rmse(t_ref, ref["v"], t_p, v_p)
    assert rmse < 1.0, f"RMSE vs MATLAB too high: {rmse:.2f}"
