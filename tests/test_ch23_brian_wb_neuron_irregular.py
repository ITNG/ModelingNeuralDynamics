from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]


def test_wb_neuron_irregular_matches_python():
    ns = load_notebook_definitions_as_module(ROOT / "brian" / "chapter23.ipynb")
    py = load_python_port(
        ROOT / "python" / "23_Entrainment_by_Excitatory_Input_Pulses" / "WB_NEURON_IRREGULAR" / "main.py"
    )

    for g_syn, t_final, v_py, _ in py.panels:
        sp_py = py.spike_times(v_py, t_final)
        sm = ns.simulate_WB_periodic_input(g_syn, t_final * ns.b2.ms)
        sp_brian = ns.spike_times_from_trace(sm.t / ns.b2.ms, sm.vm[0] / ns.b2.mV)
        # "Irregular" firing here means the trajectory sits close to a
        # bifurcation between locking regimes, so small integration
        # differences (brian's rk4 vs. this port's numpy Heun) compound
        # over the longer (800-3200ms) runs -- a looser tolerance than
        # the other ch23 tests, consistent with that sensitivity.
        assert len(sp_py) == len(sp_brian), f"g_syn={g_syn}: spike count mismatch"
        np.testing.assert_allclose(sp_brian, sp_py, atol=5.0)
