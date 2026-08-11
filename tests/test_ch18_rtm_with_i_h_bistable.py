from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "18/RTM_WITH_I_H_BISTABLE"


@pytest.mark.slow
def test_rtm_with_i_h_bistable_matches_matlab():
    # matlab overwrites v/m/h/n/r for both runs -- only the second (firing)
    # trace survives to the end of the script
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter18.ipynb")
    v_rest, v_fire, v_star, t_final, dt = ns.simulate_rtm_with_i_h_bistable()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "v_star"], timeout=60)

    assert np.isclose(v_star, ref["v_star"], atol=1e-6)
    t = np.arange(len(v_fire)) * dt
    rmse = trace_rmse(t, ref["v"], t, v_fire)
    assert rmse < 2.0, f"v RMSE vs MATLAB too high: {rmse:.4f}"
