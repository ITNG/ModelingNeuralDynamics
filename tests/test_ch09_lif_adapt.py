from pathlib import Path

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "09_Spike_Frequency_Adaptation/LIF_ADAPT"
MATLAB_DIR = "09/LIF_ADAPT"


def test_lif_adapt_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v", "w"])

    rmse_v = trace_rmse(ref["t"], ref["v"], py.t, py.v)
    rmse_w = trace_rmse(ref["t"], ref["w"], py.t, py.w)
    assert rmse_v < 0.1, f"v RMSE vs MATLAB too high: {rmse_v:.4f}"
    assert rmse_w < 0.01, f"w RMSE vs MATLAB too high: {rmse_w:.4f}"
