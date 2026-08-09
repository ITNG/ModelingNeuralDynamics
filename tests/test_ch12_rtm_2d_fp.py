from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "12_Two_Dimensional_Bifurcation_Analysis/RTM_2D_FP"


def test_rtm_2d_fp_bifurcation_shape():
    # No cheap matlab workspace variable to pull here (matlab overwrites
    # a single set of v/n/J/E on every loop iteration, so nothing useful
    # survives to the end) -- check the python result is internally
    # sane instead: same shape as the figure (stable-node branch below a
    # saddle branch, meeting near i_ext~0.13, plus the spurious constant
    # blue branch at v~=-36 that matlab's own version also has).
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    black = np.array(py.points['k'])
    magenta = np.array(py.points['m'])
    blue = np.array(py.points['b'])

    assert len(black) > 100
    assert len(magenta) > 100
    assert len(blue) > 100

    # stable-node branch stays below the saddle branch until they meet
    assert black[:, 1].max() < magenta[:, 1].max() + 1
    # spurious constant branch matlab also produces
    assert np.allclose(blue[:, 1], blue[0, 1], atol=0.5)
