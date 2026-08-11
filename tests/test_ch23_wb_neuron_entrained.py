from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "23/WB_NEURON_ENTRAINED"


@pytest.mark.slow
def test_wb_neuron_entrained_matches_matlab():
    # matlab reuses v/t_spikes for both g_syn runs -- only the second
    # (g_syn=0.14) trace survives to the end of the script
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter23.ipynb")
    v_strong, v_weak, delta = ns.simulate_wb_neuron_entrained()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "t_spikes"])

    assert np.allclose(v_weak, ref["v"], atol=1e-3)

    t_spikes_weak = ns.spike_times_from_trace(v_weak, 800.0)
    assert np.allclose(t_spikes_weak, ref["t_spikes"], atol=1e-2)
