from pathlib import Path

import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def _check(simulate_kwargs, matlab_dir, v_tol, w_tol):
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    t_p, v_p, w_p = ns.simulate_rtm_m(**simulate_kwargs)

    ref = run_matlab_script(ROOT / "matlab" / matlab_dir, "make_figure.m", ["t", "v", "w"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    rmse_w = trace_rmse(ref["t"], ref["w"], t_p, w_p)
    assert rmse_v < v_tol, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
    assert rmse_w < w_tol, f"w RMSE vs MATLAB too high: {rmse_w:.4f}"


@pytest.mark.slow
def test_rtm_m_matches_matlab():
    # v_tol looser than the resting variant: 8 full-height spikes over
    # 300ms mean even sub-ms timing jitter at the upstrokes contributes
    # non-trivial RMSE, same effect as the ch05 ERISIR trace comparison.
    _check({}, "09/RTM_M", v_tol=15.0, w_tol=0.01)


@pytest.mark.slow
def test_rtm_m_resting_matches_matlab():
    _check(dict(i_ext=0, t_final=600), "09/RTM_M_RESTING", v_tol=5.0, w_tol=0.01)
