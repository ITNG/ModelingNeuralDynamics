from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_olm_with_h_and_a_currents_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter34.ipynb")

    ref = run_matlab_script(
        ROOT / "matlab" / "34" / "OLM_WITH_H_AND_A_CURRENTS", "make_figure.m", ["t", "v", "a", "b"]
    )

    rmse_v = trace_rmse(ref["t"], ref["v"], ns.sm_ha.t / ns.b2.ms, ns.sm_ha.vm[0] / ns.b2.mV)
    rmse_ab = trace_rmse(
        ref["t"], ref["a"] * ref["b"], ns.sm_ha.t / ns.b2.ms, ns.sm_ha.a[0] * ns.sm_ha.b[0]
    )
    assert rmse_v < 1.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f} mV"
    assert rmse_ab < 0.001, f"ab RMSE vs MATLAB too high: {rmse_ab:.5f}"
