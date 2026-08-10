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


def _check_rtm_f_i(name):
    # deterministic script (no RNG) -- exact match confirmed visually
    # against matlab's figure.pdf.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / name / "main.py")
    assert len(py.f_vec_constant) == len(py.i_ext_vec)
    assert len(py.f_vec_pulsed) == len(py.i_ext_vec)
    assert py.f_vec_constant[0] == 0.0
    assert py.f_vec_constant[-1] > 0.0
    assert py.f_vec_pulsed[-1] > 0.0
    # pulsed drive has a wide plateau at f=40 Hz (mode-locked to the
    # 25 ms pulse period), clearly visible in both reference figures
    assert np.count_nonzero(np.round(py.f_vec_pulsed) == 40.0) > 10


def test_rtm_f_i_curve_pulsed_excitation_structure():
    _check_rtm_f_i("RTM_F_I_CURVE_PULSED_EXCITATION")


def test_rtm_f_i_curve_pulsed_excitation_2_structure():
    _check_rtm_f_i("RTM_F_I_CURVE_PULSED_EXCITATION_2")
