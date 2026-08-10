from pathlib import Path

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_e_to_e_heterogeneous_locks_to_common_rhythm():
    """python/24.../RTM_E_TO_E_HETEROGENEOUS draws i_ext and the pairwise
    g_syn matrix randomly, and its own comment notes this isn't meant to
    bit-match MATLAB's RNG -- only the qualitative result (a heterogeneous
    network still locks to one common oscillation) is the point. So this
    checks that outcome rather than exact spike times: most of the 30
    cells should have fired a similar number of times over the 200ms run.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter24.ipynb")
    b2 = ns.b2

    counts = [
        len(ns.spike_times_from_trace(ns.sm5.t / b2.ms, ns.sm5.vm[i] / b2.mV))
        for i in range(30)
    ]
    assert all(c > 0 for c in counts), "some cells never fired"
    assert max(counts) - min(counts) <= 2, f"spike counts too spread out to call it locked: {counts}"
