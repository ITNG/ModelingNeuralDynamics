from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_adaptation_map_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    z, phi = ns.adaptation_map()

    ref = run_matlab_script(ROOT / "matlab" / "09/ADAPTATION_MAP", "make_figure.m", ["phi"])

    # Same Heun/RK2 scheme and dt as MATLAB (this is an event/threshold-
    # crossing map, not a smooth ODE, so odeint doesn't apply) -- can use a
    # tight-ish tolerance since it's effectively the same discrete algorithm.
    assert np.allclose(phi, ref["phi"], rtol=1e-3, atol=1e-4)
