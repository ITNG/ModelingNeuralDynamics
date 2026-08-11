from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "04/HH_SOLUTION"


@pytest.mark.slow
def test_hh_solution_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter04.ipynb")
    t_p, v_p, m_p, n_p, h_p = ns.simulate_hh_solution()

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"])

    rmse = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    assert rmse < 1.0, f"RMSE vs MATLAB too high: {rmse:.2f} mV"
