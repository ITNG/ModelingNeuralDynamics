from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_plot_f_entrainment_2_fixed_point_iteration():
    # F/T/tau/epsilon match matlab's source exactly (checked by hand);
    # verify the alpha sequence is a genuine forward-iteration fixed point
    # approach (alpha_3 close to alpha_2, both inside (0, 1-epsilon))
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter23.ipynb")
    T, tau, epsilon = 25.0, 50.0, 0.6
    a1, a2, a3 = ns.simulate_f_entrainment_2(T=T, tau=tau, epsilon=epsilon)

    assert a1 == 0.0
    assert np.isclose(ns.F_entrainment(a1, T, tau, epsilon), a2)
    assert np.isclose(ns.F_entrainment(a2, T, tau, epsilon), a3)
    assert 0 < a2 < a3 < 1
