from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_rtm_m_resting_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter09.ipynb")
    from mnd.brian.core import get_step_current

    current = get_step_current(0, 600, ns.b2.ms, 0.0 * ns.b2.uA)
    sm = ns.simulate_RTM_M_neuron(current, 600 * ns.b2.ms)
    t = sm.t / ns.b2.ms
    w = sm.w[0]

    ref = run_matlab_script(ROOT / "matlab" / "09" / "RTM_M_RESTING", "make_figure.m", ["t", "w"])

    rmse = trace_rmse(ref["t"], ref["w"], t, w)
    assert rmse < 0.01, f"w RMSE vs MATLAB too high: {rmse:.4f}"
