from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_ping_5_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural/statistical properties rather than
    # exact spike times: all three panels should show substantial E and
    # I firing, and panel 1 (dense, all-to-all-ish connectivity) should
    # show a much more regular/periodic population rhythm than panel 2
    # (very sparse connectivity, near-independent firing).
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter30.ipynb")
    panel1, panel2, panel3 = ns.run_ping5_panels()
    t_e_spikes_1, i_e_spikes_1, t_i_spikes_1, i_i_spikes_1, lfp_1 = panel1
    t_e_spikes_2, i_e_spikes_2, t_i_spikes_2, i_i_spikes_2, lfp_2 = panel2
    t_e_spikes_3, i_e_spikes_3, t_i_spikes_3, i_i_spikes_3, lfp_3 = panel3

    assert len(t_e_spikes_1) > 500
    assert len(t_i_spikes_1) > 100
    assert len(t_e_spikes_2) > 500
    assert len(t_i_spikes_2) > 100
    assert len(t_e_spikes_3) > 500
    assert len(t_i_spikes_3) > 100

    t_final = 500.0

    def population_rate_cv(t_spikes, bin_width=5.0, t_final=t_final):
        bins = np.arange(0, t_final + bin_width, bin_width)
        counts, _ = np.histogram(t_spikes, bins=bins)
        return counts.std() / counts.mean()

    cv_dense = population_rate_cv(t_e_spikes_1)
    cv_sparse = population_rate_cv(t_e_spikes_2)
    assert cv_dense > cv_sparse
