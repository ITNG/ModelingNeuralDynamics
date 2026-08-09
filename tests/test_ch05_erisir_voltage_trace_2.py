import importlib.util
from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "05_The_Simple_Model_of_Neurons_in_Rodent_Brains/ERISIR_VOLTAGE_TRACE_2"
MATLAB_DIR = "05/ERISIR_VOLTAGE_TRACE_2"


def _load_main():
    path = ROOT / "python" / PYTHON_DIR / "main.py"
    spec = importlib.util.spec_from_file_location("erisir_voltage_trace_2", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_erisir_voltage_trace_2_matches_matlab():
    py = _load_main()
    t_p = py.np.arange(0, py.t_final, py.dt)
    v_p = odeint(py.derivative, py.x0, t_p)[:, 0]

    ref = run_matlab_script(
        ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"]
    )

    rmse = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    assert rmse < 5.0, f"RMSE vs MATLAB too high: {rmse:.2f} mV"
