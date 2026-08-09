from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_plot_s_matches_matlab():
    py = load_python_port(ROOT / "python" / "20_Chemical_Synapses" / "RTM_PLOT_S" / "main.py")
    py.tau_r, py.tau_d = 1.0, 2.0
    t = py.np.arange(0, py.t_final, py.dt)
    x0 = py.initial_condition(py.v)
    sol = odeint(py.derivative, x0, t)
    V, S = sol[:, 0], sol[:, -1]

    ref = run_matlab_script(ROOT / "matlab" / "20" / "RTM_PLOT_S", "make_figure.m", ["t", "v", "s"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t, V)
    rmse_s = trace_rmse(ref["t"], ref["s"], t, S)
    assert rmse_v < 2.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f} mV"
    assert rmse_s < 0.01, f"s RMSE vs MATLAB too high: {rmse_s:.4f}"
