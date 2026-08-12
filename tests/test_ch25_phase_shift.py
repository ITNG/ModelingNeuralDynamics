from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "25/PHASE_SHIFT"


@pytest.mark.slow
def test_phase_shift_matches_matlab():
    # matlab's script reuses the variable names v/t across all four
    # simulation runs, so by the end of the script they hold only the
    # last (t_pulse=120 perturbed) run's results.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter25.ipynb")
    t, v_baseline, v_perturbed = ns.simulate_phase_shift(t_pulse=120.0)
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["t", "v"], timeout=60)

    assert np.allclose(t, ref["t"], atol=1e-9)
    # spiking trace: fixed-step floating-point drift shifts spike timing
    # by a fraction of a step over 300ms, which blows up a point-by-point
    # comparison right on a spike upstroke even though the trace is
    # otherwise identical; RMSE over the whole trace is robust to that.
    rmse = trace_rmse(ref["t"], ref["v"], t, v_perturbed)
    assert rmse < 10.0
