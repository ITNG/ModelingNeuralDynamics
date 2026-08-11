from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_plot_q_smoke():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    t, V, Q, S = ns.simulate_rtm_plot_q(t_final=20.0)

    assert t.shape == V.shape == Q.shape == S.shape
    assert np.all(np.isfinite(V)) and np.all(np.isfinite(Q)) and np.all(np.isfinite(S))
    assert np.all((Q >= 0) & (Q <= 1)) and np.all((S >= 0) & (S <= 1))
