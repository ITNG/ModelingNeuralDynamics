from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]


def test_theta_firing_matches_matlab():
    """theta(t) has no reset (it's a continuous phase variable), so unlike
    QIF/LIF this one's MATLAB trace is directly comparable. make_figure.m
    doesn't keep a 't' variable (it builds the time axis inline for
    plotting), so only 'theta' is pulled from the workspace.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter08.ipynb")

    ref = run_matlab_script(ROOT / "matlab" / "08" / "THETA_FIRING", "make_figure.m", ["theta"])
    theta_ref = ref["theta"]
    t_ref = np.arange(len(theta_ref)) * 0.001

    v_ref = 1 - np.cos(theta_ref)
    v_test = 1 - np.cos(ns.theta2)

    v_test_interp = np.interp(t_ref, ns.t2, v_test)
    rmse = np.sqrt(np.mean((v_test_interp - v_ref) ** 2))
    assert rmse < 0.01, f"1-cos(theta) RMSE vs MATLAB too high: {rmse:.4f}"
