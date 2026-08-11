from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "14/ERISIR_REDUCED"


@pytest.mark.slow
def test_erisir_reduced_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter14.ipynb")
    t, v_3d, v_2d = ns.simulate_erisir_reduced()

    # matlab's script overwrites v/t -- the 2D (reduced) model's trace is
    # what's left in the workspace at the end
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"])

    rmse_2d = trace_rmse(ref["t"], ref["v"], t, v_2d)
    assert rmse_2d < 15.0, f"2D model RMSE vs MATLAB too high: {rmse_2d:.2f}"

    # sanity check for the 3D model too (no matlab reference survives for
    # it, so just check it's a plausible spike train in the same range)
    assert v_3d.max() > 30 and v_3d.min() < -70
