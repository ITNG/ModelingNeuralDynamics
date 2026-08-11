from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_notebook_definitions_as_module, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_wb_network_with_gj_subthreshold_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter21.ipynb")
    ns_python = load_notebook_definitions_as_module(ROOT / "python" / "chapter21.ipynb")

    gap_always_on = np.array([[0.0, 0.01], [0.01, 0.0]])
    v_always = ns_python.simulate_wb_network_with_gj_subthreshold(lambda v1: gap_always_on)

    def gap_subthreshold(v1, thr=-50.0):
        return gap_always_on if v1 <= thr else np.zeros((2, 2))

    v_sub = ns_python.simulate_wb_network_with_gj_subthreshold(gap_subthreshold)
    t = np.arange(v_always.shape[1]) * 0.01

    for label, v_ref, sm in [("always", v_always, ns.sm_always), ("sub", v_sub, ns.sm_sub)]:
        for cell in range(2):
            rmse = trace_rmse(t, v_ref[cell], sm.t / ns.b2.ms, sm.vm[cell] / ns.b2.mV)
            assert rmse < 1.5, f"v{cell+1} ({label}) RMSE too high: {rmse:.2f} mV"
