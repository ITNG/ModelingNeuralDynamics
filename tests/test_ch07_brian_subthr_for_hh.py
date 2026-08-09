from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_subthr_for_hh_matches_matlab():
    """Same underlying trace as LIF_NEURON_WITH_HH (see that test for the
    v(t) spike-time check). I_na/I_k inherit v's spike-upstroke RMSE
    blowup, but I_l = gl*(El-v) is smooth enough for a direct comparison.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter07.ipynb")

    ref = run_matlab_script(ROOT / "matlab" / "07" / "SUBTHR_FOR_HH", "make_figure.m", ["t", "I_l"])

    rmse = trace_rmse(ref["t"], ref["I_l"], ns.t, ns.I_l)
    assert rmse < 3.0, f"I_l RMSE vs MATLAB too high: {rmse:.2f} uA/cm^2"
