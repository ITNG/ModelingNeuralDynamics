from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_reset_threshold_matches_matlab():
    """Reuses chapter 5's WB neuron with i_ext=0.75uA -- this is the same
    trace used to pick the -67/-52mV reset/threshold pair for the LIF
    approximation used by the rest of chapter 21.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter21.ipynb")

    ref = run_matlab_script(ROOT / "matlab" / "21" / "RESET_THRESHOLD", "make_figure.m", ["t", "v"])

    rmse = trace_rmse(ref["t"], ref["v"], ns.sm_rt.t / ns.b2.ms, ns.sm_rt.vm[0] / ns.b2.mV)
    assert rmse < 1.5, f"v RMSE vs MATLAB too high: {rmse:.2f} mV"
