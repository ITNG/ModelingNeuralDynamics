from pathlib import Path

import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_lif_adapt_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    t, v, w = ns.simulate_lif_adapt()

    ref = run_matlab_script(ROOT / "matlab" / "09/LIF_ADAPT", "make_figure.m", ["t", "v", "w"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t, v)
    rmse_w = trace_rmse(ref["t"], ref["w"], t, w)
    assert rmse_v < 0.1, f"v RMSE vs MATLAB too high: {rmse_v:.4f}"
    assert rmse_w < 0.01, f"w RMSE vs MATLAB too high: {rmse_w:.4f}"
