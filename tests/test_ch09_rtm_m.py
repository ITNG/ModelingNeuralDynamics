from pathlib import Path

from scipy.integrate import odeint

from matlab_ref import load_python_port, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def _check(python_dir, matlab_dir, v_tol, w_tol):
    py = load_python_port(ROOT / "python" / python_dir / "main.py")
    t_p = py.np.arange(0, py.t_final, py.dt)
    sol = odeint(py.derivative, py.x0, t_p)
    v_p, w_p = sol[:, 0], sol[:, 3]

    ref = run_matlab_script(ROOT / "matlab" / matlab_dir, "make_figure.m", ["t", "v", "w"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    rmse_w = trace_rmse(ref["t"], ref["w"], t_p, w_p)
    assert rmse_v < v_tol, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
    assert rmse_w < w_tol, f"w RMSE vs MATLAB too high: {rmse_w:.4f}"


def test_rtm_m_matches_matlab():
    # v_tol looser than the resting variant: 8 full-height spikes over
    # 300ms mean even sub-ms timing jitter at the upstrokes contributes
    # non-trivial RMSE, same effect as the ch05 ERISIR trace comparison.
    _check(
        "09_Spike_Frequency_Adaptation/RTM_M",
        "09/RTM_M",
        v_tol=15.0, w_tol=0.01,
    )


def test_rtm_m_resting_matches_matlab():
    _check(
        "09_Spike_Frequency_Adaptation/RTM_M_RESTING",
        "09/RTM_M_RESTING",
        v_tol=5.0, w_tol=0.01,
    )
