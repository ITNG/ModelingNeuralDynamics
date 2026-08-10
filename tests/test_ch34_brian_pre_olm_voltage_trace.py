from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_pre_olm_voltage_trace_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter34.ipynb")

    ref = run_matlab_script(
        ROOT / "matlab" / "34" / "PRE_OLM_VOLTAGE_TRACE", "make_figure.m", ["t", "v"]
    )

    rmse = trace_rmse(ref["t"], ref["v"], ns.sm_pre.t / ns.b2.ms, ns.sm_pre.vm[0] / ns.b2.mV)
    assert rmse < 1.0, f"v RMSE vs MATLAB too high: {rmse:.2f} mV"
