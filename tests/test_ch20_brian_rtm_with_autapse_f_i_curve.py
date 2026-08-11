from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_brian_rtm_autapse_f_i_matches_verified_python():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    ns_python = load_notebook_definitions_as_module(ROOT / "python" / "chapter20.ipynb")
    f_forward, f_backward, I_c, I_star, i_ext_vec = ns_python.simulate_rtm_with_autapse_f_i_curve()
    result = ns.sweep_rtm_autapse(i_ext_vec)
    np.testing.assert_allclose(result["f_backward"], f_backward, atol=1.0)
    assert np.isclose(result["i_c"], I_c, atol=0.0026)
    assert np.isclose(result["i_star"], I_star, atol=0.0026)
