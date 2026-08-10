from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_notebook_as_module, load_python_port, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_wb_network_with_gj_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter21.ipynb")
    py = load_python_port(ROOT / "python" / "21_Gap_Junctions" / "WB_NETWORK_WITH_GJ" / "main.py")

    t = py.np.arange(0, py.t_final, py.dt)
    sol = odeint(py.derivative, py.x0, t)
    v1, v2 = sol[:, 0], sol[:, 1]

    rmse_v1 = trace_rmse(t, v1, ns.sm_wb.t / ns.b2.ms, ns.sm_wb.vm[0] / ns.b2.mV)
    rmse_v2 = trace_rmse(t, v2, ns.sm_wb.t / ns.b2.ms, ns.sm_wb.vm[1] / ns.b2.mV)
    assert rmse_v1 < 1.5, f"v1 RMSE vs python too high: {rmse_v1:.2f} mV"
    assert rmse_v2 < 0.05, f"v2 RMSE vs python too high: {rmse_v2:.4f} mV"
