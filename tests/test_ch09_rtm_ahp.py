from pathlib import Path

import pytest
from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def _check(python_dir, matlab_dir, v_tol, ca_tol):
    py = load_python_port(ROOT / "python" / python_dir / "main.py")
    t_p = py.np.arange(0, py.t_final, py.dt)
    sol = odeint(py.derivative, py.x0, t_p)
    v_p, ca_p = sol[:, 0], sol[:, 3]

    ref = run_matlab_script(ROOT / "matlab" / matlab_dir, "make_figure.m", ["t", "v", "ca"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    rmse_ca = trace_rmse(ref["t"], ref["ca"], t_p, ca_p)
    assert rmse_v < v_tol, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
    assert rmse_ca < ca_tol, f"ca RMSE vs MATLAB too high: {rmse_ca:.4f}"


@pytest.mark.slow
def test_rtm_ahp_matches_matlab():
    _check(
        "09_Spike_Frequency_Adaptation/RTM_AHP",
        "09/RTM_AHP",
        v_tol=15.0, ca_tol=0.01,
    )


@pytest.mark.slow
def test_rtm_ahp_resting_matches_matlab():
    _check(
        "09_Spike_Frequency_Adaptation/RTM_AHP_RESTING",
        "09/RTM_AHP_RESTING",
        v_tol=5.0, ca_tol=0.001,
    )
