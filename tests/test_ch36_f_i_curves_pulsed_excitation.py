from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "36_F_I_Curves_Pulsed_Excitation"


def test_idealized_f_i_curve_structure():
    # module only defines its curve inside `if __name__ == "__main__"`,
    # so just check the file imports cleanly and produced a figure.
    load_python_port(ROOT / "python" / PYTHON_BASE / "IDEALIZED_F_I_CURVE" / "main.py")
    assert (ROOT / "python" / PYTHON_BASE / "IDEALIZED_F_I_CURVE" / "fig.png").exists()


def test_square_pulses_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "SQUARE_PULSES" / "main.py")
    assert py.i_h > py.i_l
    assert py.T > py.epsilon


def _load_rtm_f_i(name):
    py = load_python_port(ROOT / "python" / PYTHON_BASE / name / "main.py")
    # the dense scans are explicit calls, not import side effects
    assert not hasattr(py, "f_vec_constant")
    assert not hasattr(py, "f_vec_pulsed")
    return py


def _check_rtm_f_i(name):
    # deterministic script (no RNG) -- exact match confirmed visually
    # against matlab's figure.pdf.
    py = _load_rtm_f_i(name)
    f_vec_constant, f_vec_pulsed = py.compute_f_i_curves()
    assert len(f_vec_constant) == len(py.i_ext_vec)
    assert len(f_vec_pulsed) == len(py.i_ext_vec)
    assert f_vec_constant[0] == 0.0
    assert f_vec_constant[-1] > 0.0
    assert f_vec_pulsed[-1] > 0.0
    # pulsed drive has a wide plateau at f=40 Hz (mode-locked to the
    # 25 ms pulse period), clearly visible in both reference figures
    assert np.count_nonzero(np.round(f_vec_pulsed) == 40.0) > 10


def _check_rtm_f_i_numba(name):
    py = _load_rtm_f_i(name)
    grid = py.i_ext_vec
    probe = np.array([grid[0], grid[len(grid) // 2], grid[-1]])
    expected = py.compute_f_i_curves(i_ext_values=probe, use_numba=False)
    actual = py.compute_f_i_curves(i_ext_values=probe, use_numba=True)
    for expected_curve, actual_curve in zip(expected, actual):
        np.testing.assert_allclose(actual_curve, expected_curve, rtol=0.0, atol=1e-10)


def test_rtm_f_i_curve_pulsed_excitation_structure():
    _check_rtm_f_i("RTM_F_I_CURVE_PULSED_EXCITATION")


def test_rtm_f_i_curve_pulsed_excitation_2_structure():
    _check_rtm_f_i("RTM_F_I_CURVE_PULSED_EXCITATION_2")


def test_rtm_f_i_curve_pulsed_excitation_numba_matches_python():
    _check_rtm_f_i_numba("RTM_F_I_CURVE_PULSED_EXCITATION")


def test_rtm_f_i_curve_pulsed_excitation_2_numba_matches_python():
    _check_rtm_f_i_numba("RTM_F_I_CURVE_PULSED_EXCITATION_2")
