from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "04/HH_LIMIT_CYCLE"


@pytest.mark.slow
def test_hh_limit_cycle_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter04.ipynb")
    t_p, v_p, n_p = ns.simulate_hh_limit_cycle()

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v", "n"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    rmse_n = trace_rmse(ref["t"], ref["n"], t_p, n_p)
    assert rmse_v < 1.0, f"v RMSE vs MATLAB too high: {rmse_v:.2f} mV"
    assert rmse_n < 0.01, f"n RMSE vs MATLAB too high: {rmse_n:.4f}"
