from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "04/HH_REFRACTORINESS"


@pytest.mark.slow
def test_hh_refractoriness_matches_matlab():
    """make_figure.m runs all three pulse-onset panels (-2, 5, 9) in one
    script, reusing v/t each time -- so its workspace at the end only holds
    the last (onset=9) trace. Compare that one against the matching Python
    run.
    """
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter04.ipynb")
    t_p, v_p = ns.simulate_hh_refractoriness(pulse_onset=9.0)

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"])

    rmse = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    assert rmse < 1.0, f"RMSE vs MATLAB too high: {rmse:.2f} mV"
