from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_plot_s_two_variables_smoke():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    t, V, Q, S = ns.simulate_two_stage_synapse(tau_r=10.0, tau_d=300.0, tau_d_q=10.0, t_final=50.0)

    assert t.shape == V.shape == Q.shape == S.shape
    assert np.all(np.isfinite(V)) and np.all(np.isfinite(Q)) and np.all(np.isfinite(S))
    assert np.all((Q >= 0) & (Q <= 1)) and np.all((S >= 0) & (S <= 1))
