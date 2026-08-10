from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_square_waves_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter19.ipynb")
    sm = ns.simulate_inapik_slow_k(g_k_slow=5.0 * ns.b2.msiemens)
    t = sm.t / ns.b2.ms
    v = sm.v[0] / ns.b2.mV

    ref = run_matlab_script(
        ROOT / "matlab" / "19" / "SQUARE_WAVES", "make_figure.m", ["t", "v"]
    )

    assert trace_rmse(ref["t"], ref["v"], t, v) < 0.01


def test_square_waves_highlights_the_slow_quasi_steady_segments():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter19.ipynb")
    sm = ns.simulate_inapik_slow_k(g_k_slow=5.0 * ns.b2.msiemens)
    t = sm.t / ns.b2.ms
    v = sm.v[0] / ns.b2.mV
    dt = float(t[1] - t[0])

    v_left, v_right = v[:-2], v[2:]
    ind = np.where(np.abs(v_right - v_left) / dt < 1)[0]

    assert 0 < len(ind) < len(v)
    # the highlighted samples sit on the slow inter-burst ramp, well below
    # the fast spike upstrokes
    assert np.all(v[ind + 1] < -30)
