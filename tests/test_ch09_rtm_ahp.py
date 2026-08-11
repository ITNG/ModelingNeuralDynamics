from pathlib import Path

import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def _check(simulate_kwargs, matlab_dir, v_tol, ca_tol):
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    t_p, v_p, ca_p = ns.simulate_rtm_ahp(**simulate_kwargs)

    ref = run_matlab_script(ROOT / "matlab" / matlab_dir, "make_figure.m", ["t", "v", "ca"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    rmse_ca = trace_rmse(ref["t"], ref["ca"], t_p, ca_p)
    assert rmse_v < v_tol, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
    assert rmse_ca < ca_tol, f"ca RMSE vs MATLAB too high: {rmse_ca:.4f}"


@pytest.mark.slow
def test_rtm_ahp_matches_matlab():
    _check({}, "09/RTM_AHP", v_tol=15.0, ca_tol=0.01)


@pytest.mark.slow
def test_rtm_ahp_resting_matches_matlab():
    _check(dict(i_ext=0, t_final=600), "09/RTM_AHP_RESTING", v_tol=5.0, ca_tol=0.001)
