from pathlib import Path

from matlab_ref import load_notebook_as_module, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_v_v_tilde_matches_matlab():
    """make_figure.m reuses its 'v' array for both traces (the tilde-w_k
    run overwrites the plain run in the workspace), so this compares
    against the MATLAB update rule reproduced directly for each w_k
    rather than pulling both out of one script run.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter09.ipynb")

    tau_m, I, tau_w = 10.0, 0.13, 40.0
    dt = 0.01

    def matlab_v(wk):
        import numpy as np

        m_steps = round(100 / dt)
        v = np.zeros(m_steps + 1)
        for k in range(m_steps):
            t = k * dt
            v_inc = -v[k] / tau_m + I - wk * np.exp(-t / tau_w) * v[k]
            v_tmp = v[k] + dt / 2 * v_inc
            v_inc = -v_tmp / tau_m + I - wk * np.exp(-(t + dt / 2) / tau_w) * v_tmp
            v[k + 1] = v[k] + dt * v_inc
        t_arr = np.arange(m_steps + 1) * dt
        return t_arr, v

    t_ref, v_ref = matlab_v(0.05)
    t_ref_tilde, v_tilde_ref = matlab_v(0.08)

    rmse_v = trace_rmse(t_ref, v_ref, ns.sm10.t / ns.b2.ms, ns.sm10.v[0])
    rmse_vt = trace_rmse(t_ref_tilde, v_tilde_ref, ns.sm10_tilde.t / ns.b2.ms, ns.sm10_tilde.v[0])

    assert rmse_v < 0.001, f"v RMSE vs MATLAB too high: {rmse_v:.5f}"
    assert rmse_vt < 0.001, f"v_tilde RMSE vs MATLAB too high: {rmse_vt:.5f}"
