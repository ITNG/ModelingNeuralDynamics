from pathlib import Path

from matlab_ref import load_notebook_as_module, load_notebook_definitions_as_module, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_wb_network_with_gj_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter21.ipynb")
    ns_python = load_notebook_definitions_as_module(ROOT / "python" / "chapter21.ipynb")

    t, v1, v2 = ns_python.simulate_wb_network_with_gj()

    rmse_v1 = trace_rmse(t, v1, ns.sm_wb.t / ns.b2.ms, ns.sm_wb.vm[0] / ns.b2.mV)
    rmse_v2 = trace_rmse(t, v2, ns.sm_wb.t / ns.b2.ms, ns.sm_wb.vm[1] / ns.b2.mV)
    assert rmse_v1 < 1.5, f"v1 RMSE vs python too high: {rmse_v1:.2f} mV"
    assert rmse_v2 < 0.05, f"v2 RMSE vs python too high: {rmse_v2:.4f} mV"
