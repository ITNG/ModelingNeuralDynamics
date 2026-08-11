from pathlib import Path

import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_v_v_tilde_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    t, v, v_tilde = ns.v_v_tilde()

    # matlab's make_figure.m overwrites its own `v` in place: first with
    # w_k, then again with tilde_w_k. Only the second (tilde_v) survives to
    # be saved -- so that's the only trace matlab and python have both in
    # scope with matching semantics.
    ref = run_matlab_script(ROOT / "matlab" / "09/V_V_TILDE", "make_figure.m", ["v"])

    rmse = trace_rmse(t, ref["v"], t, v_tilde)
    assert rmse < 0.05, f"RMSE vs MATLAB too high: {rmse:.4f}"
