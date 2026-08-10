from pathlib import Path

from matlab_ref import load_notebook_as_module, load_python_port, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_wb_network_with_gj_subthreshold_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter21.ipynb")
    py = load_python_port(
        ROOT / "python" / "21_Gap_Junctions" / "WB_NETWORK_WITH_GJ_SUBTHRESHOLD" / "main.py"
    )
    t = py.np.arange(py.m_steps + 1) * py.dt

    for label, v_ref, sm in [("always", py.v_always, ns.sm_always), ("sub", py.v_sub, ns.sm_sub)]:
        for cell in range(2):
            rmse = trace_rmse(t, v_ref[cell], sm.t / ns.b2.ms, sm.vm[cell] / ns.b2.mV)
            assert rmse < 1.5, f"v{cell+1} ({label}) RMSE too high: {rmse:.2f} mV"
