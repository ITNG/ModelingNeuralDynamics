from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_ellipses_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter19.ipynb")
    sm = ns.simulate_erisir_slow_k(g_k_slow=1.5 * ns.b2.msiemens)
    t = sm.t / ns.b2.ms
    v = sm.v[0] / ns.b2.mV

    ref = run_matlab_script(
        ROOT / "matlab" / "19" / "ELLIPSES", "make_figure.m", ["t", "v"]
    )

    assert trace_rmse(ref["t"], ref["v"], t, v) < 0.01


def test_ellipses_highlights_local_maxima_and_minima_envelope():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter19.ipynb")
    sm = ns.simulate_erisir_slow_k(g_k_slow=1.5 * ns.b2.msiemens)
    v = sm.v[0] / ns.b2.mV

    v_left, v_center, v_right = v[:-2], v[1:-1], v[2:]
    maxima = np.where((v_center > v_left) & (v_center > v_right))[0]
    minima = np.where((v_center < v_left) & (v_center < v_right))[0]

    assert len(maxima) > 0
    assert len(minima) > 0
    # the envelope spans from the spike troughs up to the spike peaks (small
    # ripples during the slow interburst phase also count as local extrema,
    # same as the original MATLAB/Python highlighting logic)
    assert v[maxima + 1].max() > 40
    assert v[minima + 1].min() < -80
