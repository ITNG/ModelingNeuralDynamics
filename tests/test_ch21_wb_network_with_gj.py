from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_wb_network_with_gj_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter21.ipynb")
    t, v1, v2 = ns.simulate_wb_network_with_gj()

    ref = run_matlab_script(ROOT / "matlab" / "21" / "WB_NETWORK_WITH_GJ", "make_figure.m", ["t", "v"])

    rmse_v1 = trace_rmse(ref["t"], ref["v"][0], t, v1)
    rmse_v2 = trace_rmse(ref["t"], ref["v"][1], t, v2)
    assert rmse_v1 < 1.5, f"v1 RMSE vs MATLAB too high: {rmse_v1:.2f} mV"
    assert rmse_v2 < 0.05, f"v2 RMSE vs MATLAB too high: {rmse_v2:.4f} mV"
