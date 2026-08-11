from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_e_to_e_network_smoke():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter24.ipynb")
    t_spikes, i_spikes = ns.simulate_rtm_e_to_e_network(t_final=50.0)

    assert len(t_spikes) == len(i_spikes)
    assert len(t_spikes) > 0
    assert np.all((i_spikes >= 1) & (i_spikes <= 30))
