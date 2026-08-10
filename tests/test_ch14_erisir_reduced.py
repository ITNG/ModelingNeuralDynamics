from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "14_Model_Neurons_of_Bifurcation_Type_2/ERISIR_REDUCED"
MATLAB_DIR = "14/ERISIR_REDUCED"


@pytest.mark.slow
def test_erisir_reduced_matches_matlab():
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    t = py.np.arange(0, py.t_final, py.dt)
    v_3d = odeint(py.derivative_3d, py.x0_3d, t)[:, 0]
    v_2d = odeint(py.derivative_2d, py.x0_2d, t)[:, 0]

    # matlab's script overwrites v/t -- the 2D (reduced) model's trace is
    # what's left in the workspace at the end
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["t", "v"])

    rmse_2d = trace_rmse(ref["t"], ref["v"], t, v_2d)
    assert rmse_2d < 15.0, f"2D model RMSE vs MATLAB too high: {rmse_2d:.2f}"

    # sanity check for the 3D model too (no matlab reference survives for
    # it, so just check it's a plausible spike train in the same range)
    assert v_3d.max() > 30 and v_3d.min() < -70
