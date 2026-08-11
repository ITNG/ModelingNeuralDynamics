from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_erisir_2d_fp_bifurcation_shape():
    # Same situation as ch12/RTM_2D_FP: matlab overwrites a single v/n/J/E
    # on every loop iteration, so nothing useful survives to compare
    # numerically. Check the python result is internally sane instead --
    # a fold with a stable branch (black/red) below an unstable saddle
    # branch (magenta), meeting near i_ext~6.2, plus the spurious constant
    # unstable-node branch at v~=-14 that matlab's own version also has.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter14.ipynb")
    points, i_ext_vec = ns.simulate_erisir_2d_fp()

    black = np.array(points['k'])
    red = np.array(points['r'])
    magenta = np.array(points['m'])
    blue = np.array(points['b'])

    assert len(black) > 100
    assert len(red) > 100
    assert len(magenta) > 100
    assert len(blue) > 100

    # stable branch (node then focus) stays below the saddle branch
    assert max(black[:, 1].max(), red[:, 1].max()) < magenta[:, 1].max()
    # spurious constant branch matlab also produces
    assert np.allclose(blue[:, 1], blue[0, 1], atol=0.5)
