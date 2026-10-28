from pathlib import Path

import numpy as np
import pytest

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


@pytest.mark.parametrize("name,kwargs", [
    ("spiral_sub", dict(I=.1, r0=.3, theta0=0., t_final=5., dt=.01)),
    ("spiral_sub_2", dict(I=-.2, r0=1.1, theta0=0., t_final=5., dt=.001, reverse=True)),
    ("spiral_sup", dict(I=.5, r0=1.2, theta0=0., t_final=3., dt=.0001, reverse=True)),
])
def test_escaping_spirals_stop_at_plot_domain_without_changing_inner_path(name, kwargs):
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter13.ipynb")
    spiral = getattr(ns, name)
    with np.errstate(over="raise", invalid="raise"):
        x, y = spiral(**kwargs, max_radius=2.)
        assert np.isfinite(x).all() and np.isfinite(y).all()
        assert len(x) < round(kwargs["t_final"]/kwargs["dt"])+1
        short = {**kwargs, "t_final": .02}
        np.testing.assert_array_equal(spiral(**short, max_radius=2.), spiral(**short))
