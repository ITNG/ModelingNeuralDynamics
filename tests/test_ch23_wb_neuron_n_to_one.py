from pathlib import Path

import numpy as np

from matlab_ref import load_python_port, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "23_Entrainment_by_Excitatory_Input_Pulses/WB_NEURON_N_TO_ONE"
MATLAB_DIR = "23/WB_NEURON_N_TO_ONE"


@pytest.mark.slow
def test_wb_neuron_n_to_one_matches_matlab():
    # matlab reuses v/t_spikes for both g_syn runs -- only the second
    # (g_syn=0.15, t_final=1600) trace survives to the end of the script
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v", "t_spikes"],
                             timeout=60)

    g_syn, t_final, v, delta = py.panels[-1]
    assert g_syn == 0.15 and t_final == 1600.
    assert np.allclose(v, ref["v"], atol=1e-3)

    t_spikes = py.spike_times(v, t_final)
    assert np.allclose(t_spikes, ref["t_spikes"], atol=1e-2)
