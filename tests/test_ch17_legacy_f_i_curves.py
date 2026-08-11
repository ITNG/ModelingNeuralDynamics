from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]

# The three legacy full-model examples (HH_F_I_CURVE, RTM_F_I_CURVE,
# RTM_F_I_CURVE_AT_ONSET) never had MATLAB-comparison test coverage even
# before this notebook existed -- kept as reference material, not verified
# results. These smoke tests use small scans to stay fast.


def test_chapter17_notebook_hh_f_i_curve_legacy():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    f_forward, f_backward, i_ext_vec = ns.simulate_hh_f_i_curve_legacy(
        i_ext_vec=np.array([3.0, 8.0, 13.0]), t_final=100.0
    )
    assert f_forward.shape == f_backward.shape == i_ext_vec.shape == (3,)
    assert np.all(np.isfinite(f_forward)) and np.all(np.isfinite(f_backward))
    assert f_forward[0] == f_backward[0] == 0.0
    assert f_forward[-1] > 70.0 and f_backward[-1] > 70.0


def test_chapter17_notebook_rtm_f_i_curve_legacy():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    f_forward, f_backward, i_ext_vec = ns.simulate_rtm_f_i_curve_legacy(
        i_ext_vec=np.array([0.0, 0.5, 1.0]), t_final=100.0
    )
    assert f_forward.shape == f_backward.shape == i_ext_vec.shape == (3,)
    assert np.all(np.isfinite(f_forward)) and np.all(np.isfinite(f_backward))
    assert f_forward[0] == f_backward[0] == 0.0
    assert f_forward[-1] > 40.0 and f_backward[-1] > 40.0


def test_chapter17_notebook_rtm_f_i_curve_at_onset_legacy():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter17.ipynb")
    f_forward, i_ext_vec, I_c, C = ns.simulate_rtm_f_i_curve_at_onset_legacy(
        n_points=3, t_final=100.0
    )
    assert f_forward.shape == i_ext_vec.shape == (3,)
    assert np.array_equal(f_forward, np.zeros(3))
    assert I_c is None and C is None
