from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_e_to_e_heterogeneous_is_sane():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter24.ipynb")
    t_spikes, i_spikes, i_ext, g_syn = ns.simulate_rtm_e_to_e_heterogeneous(
        N=4, t_final=100.0, seed=63806
    )
    N = len(i_ext)

    assert N == 4
    assert i_ext.shape == (N,)
    assert np.all((0.25 <= i_ext) & (i_ext <= 0.35))
    assert g_syn.shape == (N, N)
    assert np.all(np.diag(g_syn) == 0)
    off_diag = g_syn[~np.eye(N, dtype=bool)]
    assert np.all((0.00625 <= off_diag) & (off_diag <= 0.00875))
    counts = np.array([np.sum(i_spikes == i) for i in range(1, N + 1)])
    assert np.array_equal(counts, np.array([2, 2, 2, 2]))
