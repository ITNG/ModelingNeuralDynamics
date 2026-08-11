from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "05/ERISIR_VOLTAGE_TRACE_2"


@pytest.mark.slow
def test_erisir_voltage_trace_2_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter05.ipynb")
    t_p, v_p = ns.simulate_erisir_voltage_trace(n_power=4)

    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"]
    )

    rmse = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    assert rmse < 5.0, f"RMSE vs MATLAB too high: {rmse:.2f} mV"
