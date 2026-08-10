from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_tau_m_for_hh_matches_matlab():
    """tau(t) is a smooth function of m/n/h -- no spike upstrokes -- so a
    plain whole-trace RMSE is reliable here, unlike v(t) itself.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter07.ipynb")

    ref = run_matlab_script(ROOT / "matlab" / "07" / "TAU_M_FOR_HH", "make_figure.m", ["t", "tau"])

    rmse = trace_rmse(ref["t"], ref["tau"], ns.t, ns.tau)
    assert rmse < 0.5, f"tau RMSE vs MATLAB too high: {rmse:.3f} ms"
