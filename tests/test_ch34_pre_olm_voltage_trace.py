from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "34/PRE_OLM_VOLTAGE_TRACE"


@pytest.mark.slow
def test_pre_olm_voltage_trace_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter34.ipynb")
    t, v = ns.simulate_pre_olm_voltage_trace()

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"])

    rmse = trace_rmse(ref["t"], ref["v"], t, v)
    assert rmse < 1.0, f"RMSE vs MATLAB too high: {rmse:.2f} mV"
