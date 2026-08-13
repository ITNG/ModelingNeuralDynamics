from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "29/RIVER"


@pytest.mark.slow
def test_river_matches_matlab():
    # matlab's script reuses theta/g_syn for many trajectories, so by the
    # end of the script they hold only the last-computed trajectory
    # (theta(1)=1.25, g_syn(1)=0.25, integrated backwards).
    #
    # NOTE: this trajectory starts near the separatrix/unstable manifold
    # and is integrated backward for ~6200 steps -- tiny floating-point
    # differences between MATLAB's and Python's Heun integrators compound
    # over that many steps and the late trajectory can diverge to a
    # different basin (confirmed: this same divergence occurs against the
    # original pre-notebook main.py script too, so it is not a porting
    # regression). g_ast (a stable/well-conditioned quantity) matches
    # MATLAB tightly; theta_unstable/g_syn_unstable's tight atol=1e-6
    # full-trajectory comparison is inherently fragile for a
    # separatrix-adjacent backward integration and may legitimately fail.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter29.ipynb")
    river = ns.simulate_river()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["theta", "g_syn", "g_ast"], timeout=30)

    assert np.isclose(river['g_ast'], ref["g_ast"], atol=1e-3)
    assert np.allclose(river['theta_unstable'], ref["theta"], atol=1e-6)
    assert np.allclose(river['g_syn_unstable'], ref["g_syn"], atol=1e-6)
