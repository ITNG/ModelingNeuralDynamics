from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "23_Entrainment_by_Excitatory_Input_Pulses/PLOT_F_ENTRAINMENT_2"


def test_plot_f_entrainment_2_fixed_point_iteration():
    # F/T/tau/epsilon match matlab's source exactly (checked by hand);
    # verify the alpha sequence is a genuine forward-iteration fixed point
    # approach (alpha_3 close to alpha_2, both inside (0, 1-epsilon))
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    a1, a2, a3 = py.alpha
    assert a1 == 0.
    assert np.isclose(py.F(a1), a2)
    assert np.isclose(py.F(a2), a3)
    assert 0 < a2 < a3 < 1
