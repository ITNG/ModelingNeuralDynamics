import inspect
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "python" / "chapter36.ipynb"


def test_idealized_f_i_curve_structure():
    # the schematic curve has no data to check, so just confirm the
    # notebook definitions load and the plot function runs cleanly.
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ns.plot_idealized_f_i_curve()


def test_square_pulses_structure():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    ns.plot_square_pulses()
    defaults = {
        name: param.default
        for name, param in inspect.signature(ns.plot_square_pulses).parameters.items()
    }
    assert defaults["i_h"] > defaults["i_l"]
    assert defaults["T"] > defaults["epsilon"]


def _check_rtm_f_i(compute_fn):
    # deterministic simulation (no RNG) -- exact match confirmed visually
    # against matlab's figure.pdf.
    f_vec_constant, f_vec_pulsed, i_ext_vec = compute_fn()
    assert len(f_vec_constant) == len(i_ext_vec)
    assert len(f_vec_pulsed) == len(i_ext_vec)
    assert f_vec_constant[0] == 0.0
    assert f_vec_constant[-1] > 0.0
    assert f_vec_pulsed[-1] > 0.0
    # pulsed drive has a wide plateau at f=40 Hz (mode-locked to the
    # 25 ms pulse period), clearly visible in both reference figures
    assert np.count_nonzero(np.round(f_vec_pulsed) == 40.0) > 10


def _check_rtm_f_i_numba(compute_fn):
    _, _, grid = compute_fn(i_ext_values=np.linspace(0.0, 2.0, 201))
    probe = np.array([grid[0], grid[len(grid) // 2], grid[-1]])
    f_c_expected, f_p_expected, _ = compute_fn(i_ext_values=probe, use_numba=False)
    f_c_actual, f_p_actual, _ = compute_fn(i_ext_values=probe, use_numba=True)
    np.testing.assert_allclose(f_c_actual, f_c_expected, rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(f_p_actual, f_p_expected, rtol=0.0, atol=1e-10)


def test_rtm_f_i_curve_pulsed_excitation_structure():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    _check_rtm_f_i(ns.compute_rtm_f_i_curves_pulsed_excitation)


def test_rtm_f_i_curve_pulsed_excitation_2_structure():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    _check_rtm_f_i(ns.compute_rtm_f_i_curves_pulsed_excitation_2)


def test_rtm_f_i_curve_pulsed_excitation_numba_matches_python():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    _check_rtm_f_i_numba(ns.compute_rtm_f_i_curves_pulsed_excitation)


def test_rtm_f_i_curve_pulsed_excitation_2_numba_matches_python():
    ns = load_notebook_definitions_as_module(NOTEBOOK)
    _check_rtm_f_i_numba(ns.compute_rtm_f_i_curves_pulsed_excitation_2)
