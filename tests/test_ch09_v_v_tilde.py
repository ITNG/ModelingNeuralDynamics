from pathlib import Path

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "09_Spike_Frequency_Adaptation/V_V_TILDE"
MATLAB_DIR = "09/V_V_TILDE"


def test_v_v_tilde_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    # matlab's make_figure.m overwrites its own `v` in place: first with
    # w_k, then again with tilde_w_k. Only the second (tilde_v) survives to
    # be saved -- so that's the only trace matlab and python have both in
    # scope with matching semantics.
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v"])

    rmse = trace_rmse(py.t, ref["v"], py.t, py.v_tilde)
    assert rmse < 0.05, f"RMSE vs MATLAB too high: {rmse:.4f}"
