from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_wb_neuron_entrained_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter23.ipynb")
    ns_python = load_notebook_definitions_as_module(ROOT / "python" / "chapter23.ipynb")
    v_strong, v_weak, delta = ns_python.simulate_wb_neuron_entrained()

    for g_syn, sm, v_py in [(0.195, ns.sm_strong, v_strong), (0.14, ns.sm_weak, v_weak)]:
        sp_py = ns_python.spike_times_from_trace(v_py, 800.0)
        sp_brian = ns.spike_times_from_trace(sm.t / ns.b2.ms, sm.vm[0] / ns.b2.mV)
        assert len(sp_py) == len(sp_brian), f"g_syn={g_syn}: spike count mismatch"
        np.testing.assert_allclose(sp_brian, sp_py, atol=1.0)
