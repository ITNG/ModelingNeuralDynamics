from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "24/RTM_E_TO_E_NETWORK_2"


@pytest.mark.slow
def test_rtm_e_to_e_network_2_matches_matlab():
    # numba-jitted: the uncompiled 10^6-step, N=30 sweep took minutes;
    # this now finishes in seconds.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter24.ipynb")
    t_spikes, i_spikes = ns.simulate_rtm_e_to_e_network_2()
    ind1 = np.where(i_spikes == 1)[0]
    frequency = 1000 / (t_spikes[ind1[-1]] - t_spikes[ind1[-2]])

    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["frequency"], timeout=300)

    assert np.isclose(frequency, ref["frequency"], atol=0.5)
