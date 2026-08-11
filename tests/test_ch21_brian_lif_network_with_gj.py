from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_notebook_definitions_as_module, spike_times

ROOT = Path(__file__).resolve().parents[1]


def test_lif_network_with_gj_matches_python():
    """python/chapter21.ipynb is already MATLAB-verified, so compare
    against it directly rather than re-running MATLAB -- much cheaper,
    and it's the same numerical target.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter21.ipynb")
    ns_python = load_notebook_definitions_as_module(ROOT / "python" / "chapter21.ipynb")

    v_c, spikes_c = ns_python.simulate_lif_network_with_gj(3.0 * 0.01)
    v_u, spikes_u = ns_python.simulate_lif_network_with_gj(0.0)

    for k in range(2):
        b_c = spike_times(ns.sm_c.t / ns.b2.ms, ns.sm_c.v[k], threshold=0.999)
        b_u = spike_times(ns.sm_u.t / ns.b2.ms, ns.sm_u.v[k], threshold=0.999)
        assert len(b_c) == len(spikes_c[k])
        assert len(b_u) == len(spikes_u[k])
        np.testing.assert_allclose(b_c, spikes_c[k], atol=0.1)
        np.testing.assert_allclose(b_u, spikes_u[k], atol=0.1)
