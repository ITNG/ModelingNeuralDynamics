from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "15_Canard_Explosions/CANARD"
MATLAB_DIR = "15/CANARD"


@pytest.mark.slow
def test_canard_matches_matlab():
    # ~70s in matlab, ~75s in python (both bisect on I for 5 target
    # amplitudes, each requiring a 1e6-step trajectory) -- expensive but
    # deliberately so, this is the point of the "canard explosion" figure.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["I", "amp"], timeout=180)

    # matlab reuses I/amp across the ijk loop -- only the last (amp_target
    # =3.78) survives to the end of the script
    amp_target, I, v, n = py.results[-1]
    assert amp_target == 3.78
    assert np.isclose(ref["amp"], 3.78, atol=1e-6)
    assert np.isclose(I, ref["I"], atol=1e-4)
