from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]


def test_brian_rtm_autapse_f_i_matches_verified_python():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter20.ipynb")
    py = load_python_port(
        ROOT
        / "python"
        / "20_Chemical_Synapses"
        / "RTM_WITH_AUTAPSE_F_I_CURVE"
        / "main.py"
    )
    result = ns.sweep_rtm_autapse(py.i_ext_vec)
    np.testing.assert_allclose(result["f_backward"], py.f_backward, atol=1.0)
    assert np.isclose(result["i_c"], py.I_c, atol=0.0026)
    assert np.isclose(result["i_star"], py.I_star, atol=0.0026)
