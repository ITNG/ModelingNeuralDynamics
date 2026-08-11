from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_chapter13_notebook_defines_working_examples():
    """Fast (non-MATLAB) sanity check that runs on every push/PR.

    Chapter 13 is all closed-form radial-flow schematics; test_ch13_hopf_sub.py
    is the only MATLAB-comparison test and doesn't touch the python module at
    all (it recomputes the closed-form r vector inline), so it needed no
    repointing. This just exercises the three normal-form spiral integrators
    for shape/finiteness.
    """
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter13.ipynb")

    for spiral, kwargs in [
        (ns.spiral_sub, dict(I=-0.1, r0=0.1, theta0=0.0, t_final=5)),
        (ns.spiral_sub_2, dict(I=-0.2, r0=0.4, theta0=0.0, t_final=5)),
        (ns.spiral_sup, dict(I=0.5, r0=0.4, theta0=0.0, t_final=3, dt=0.001)),
    ]:
        x, y = spiral(**kwargs)
        assert x.shape == y.shape
        assert np.all(np.isfinite(x)) and np.all(np.isfinite(y))
