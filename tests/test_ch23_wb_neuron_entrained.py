from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "23_Entrainment_by_Excitatory_Input_Pulses/WB_NEURON_ENTRAINED"
MATLAB_DIR = "23/WB_NEURON_ENTRAINED"


def test_wb_neuron_entrained_matches_matlab():
    # matlab reuses v/t_spikes for both g_syn runs -- only the second
    # (g_syn=0.14) trace survives to the end of the script
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "t_spikes"])

    assert np.allclose(py.v_weak, ref["v"], atol=1e-3)

    t_spikes_weak = py.spike_times(py.v_weak)
    assert np.allclose(t_spikes_weak, ref["t_spikes"], atol=1e-2)
