from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_s_buildup_smoke():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    t, V, S, period = ns.simulate_s_buildup(tau_r=10.0, tau_d=300.0, tau_d_q=5.0, i_ext=0.2)

    assert np.all(np.isfinite(V)) and np.all(np.isfinite(S))
    assert period > 0
