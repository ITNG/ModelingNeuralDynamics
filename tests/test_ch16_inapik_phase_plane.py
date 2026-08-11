from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "16/INAPIK_PHASE_PLANE"


@pytest.mark.slow
def test_inapik_phase_plane_matches_matlab():
    # matlab overwrites v/n/t_final/dt on every panel -- only the last
    # (I=4.4) trajectory survives to the end of the script
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter16.ipynb")
    panels = ns.simulate_inapik_phase_plane()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "n"])

    p = panels[-1]
    t_ref = np.arange(len(ref["v"])) * p['dt']
    t_py = np.arange(len(p['v_full'])) * p['dt']
    rmse_v = trace_rmse(t_ref, ref["v"], t_py, p['v_full'])
    rmse_n = trace_rmse(t_ref, ref["n"], t_py, p['n_full'])
    assert rmse_v < 5.0, f"v RMSE vs MATLAB too high: {rmse_v:.3f}"
    assert rmse_n < 0.05, f"n RMSE vs MATLAB too high: {rmse_n:.4f}"
