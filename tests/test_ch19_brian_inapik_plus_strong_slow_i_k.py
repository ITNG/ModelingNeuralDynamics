from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_inapik_plus_strong_slow_i_k_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter19.ipynb")
    sm = ns.simulate_inapik_slow_k(g_k_slow=20.0 * ns.b2.msiemens)
    t = sm.t / ns.b2.ms
    v = sm.v[0] / ns.b2.mV

    ref = run_matlab_script(
        ROOT / "matlab" / "19" / "INAPIK_PLUS_STRONG_SLOW_I_K",
        "make_figure.m", ["t", "v"],
    )

    assert trace_rmse(ref["t"], ref["v"], t, v) < 0.01
