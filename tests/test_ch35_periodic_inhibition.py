from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "35_Periodic_Inhibition"


def test_oscillations_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "OSCILLATIONS" / "main.py")
    phi0 = py.phi_of(1e-5)
    phi10 = py.phi_of(10)
    # larger alpha sharpens the pulses: same mean (normalized to 1 in
    # the sense that mean over one period is fixed), but higher peak
    assert phi10.max() > phi0.max()


def test_periodic_inhibition_structure():
    # deterministic script (no RNG) -- exact match confirmed visually
    # against matlab's figure.pdf, so check the exact spike frequencies.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PERIODIC_INHIBITION" / "main.py")
    _, _, freq_015 = py.run(0.15, use_periodic=True)
    _, _, freq_02 = py.run(0.2, use_periodic=True)
    assert freq_015 == 40.0
    assert freq_02 == 80.0
    # tonic (mean) inhibition never lets I=0.15 or 0.2 reach threshold
    _, _, freq_bar_015 = py.run(0.15, use_periodic=False)
    assert freq_bar_015 == 0.0


def test_periodic_inhibition_3_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PERIODIC_INHIBITION_3" / "main.py")
    _, _, freq_015 = py.run(0.15, use_periodic=True)
    _, _, freq_02 = py.run(0.2, use_periodic=True)
    assert freq_015 == 40.0
    assert freq_02 == 80.0


def test_periodic_inhibition_2_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties instead of exact spikes.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PERIODIC_INHIBITION_2" / "main.py")
    _, _, freq = py.run(use_periodic=True)
    _, _, freq_bar = py.run(use_periodic=False)
    assert freq > freq_bar


def test_periodic_inhibition_f_i_curve_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PERIODIC_INHIBITION_F_I_CURVE" / "main.py")
    # the dense scan is an explicit call, not an import side effect
    assert not hasattr(py, "f_vec")
    f_vec = py.compute_f_i_curve()
    # visually confirmed step values (40, 80, 120, 160 Hz plateaus)
    assert set(np.round(np.unique(f_vec))) >= {0.0, 40.0, 80.0, 120.0, 160.0}
    assert f_vec[0] == 0.0
    assert f_vec[-1] == 160.0


def test_periodic_inhibition_f_i_curve_numba_matches_python():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PERIODIC_INHIBITION_F_I_CURVE" / "main.py")
    grid = py.I_vec
    probe = np.array([grid[0], grid[len(grid) // 2], grid[-1]])
    expected = py.compute_f_i_curve(i_values=probe, use_numba=False)
    actual = py.compute_f_i_curve(i_values=probe, use_numba=True)
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-10)


def test_periodic_inhibition_f_i_curve_2_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PERIODIC_INHIBITION_F_I_CURVE_2" / "main.py")
    assert set(np.round(np.unique(py.f_vec))) >= {0.0, 40.0, 80.0, 120.0, 160.0}


def _load_rtm_inhibition(name):
    py = load_python_port(ROOT / "python" / PYTHON_BASE / name / "main.py")
    # the dense scans are explicit calls, not import side effects
    assert not hasattr(py, "f_vec_tonic")
    assert not hasattr(py, "f_vec_periodic")
    return py


def _check_rtm_inhibition_numba(name):
    py = _load_rtm_inhibition(name)
    grid = py.i_ext_vec
    probe = np.array([grid[0], grid[len(grid) // 2], grid[-1]])
    expected = py.compute_f_i_curves(i_ext_values=probe, use_numba=False)
    actual = py.compute_f_i_curves(i_ext_values=probe, use_numba=True)
    for expected_curve, actual_curve in zip(expected, actual):
        np.testing.assert_allclose(actual_curve, expected_curve, rtol=0.0, atol=1e-10)


def test_rtm_f_i_curve_with_inhibition_structure():
    py = _load_rtm_inhibition("RTM_F_I_CURVE_WITH_INHIBITION")
    f_vec_tonic, f_vec_periodic = py.compute_f_i_curves()
    assert len(f_vec_tonic) == len(py.i_ext_vec)
    assert len(f_vec_periodic) == len(py.i_ext_vec)
    # both curves are silent at low drive and fire at high drive
    assert f_vec_tonic[0] == 0.0
    assert f_vec_tonic[-1] > 0.0
    assert f_vec_periodic[-1] > 0.0


def test_rtm_f_i_curve_with_inhibition_numba_matches_python():
    _check_rtm_inhibition_numba("RTM_F_I_CURVE_WITH_INHIBITION")


def test_rtm_f_i_curve_with_inhibition_2_structure():
    py = _load_rtm_inhibition("RTM_F_I_CURVE_WITH_INHIBITION_2")
    _, f_vec_periodic = py.compute_f_i_curves()
    assert len(f_vec_periodic) == len(py.i_ext_vec)
    assert f_vec_periodic[-1] > 0.0


def test_rtm_f_i_curve_with_inhibition_2_numba_matches_python():
    _check_rtm_inhibition_numba("RTM_F_I_CURVE_WITH_INHIBITION_2")
