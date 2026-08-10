from pathlib import Path

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_two_cell_network_locks_to_nonzero_offset():
    """python/24.../RTM_TWO_CELL_NETWORK shows that a large-enough
    perturbation doesn't decay back to synchrony -- the pair settles into
    a persistent, non-zero phase-locked offset around +-1.5ms (confirmed
    against the python port's own t2-t1 array, which itself settles near
    +1.5ms). Brian's version settles near -1.5ms: in a symmetric
    two-oscillator system, which of the two symmetric locked branches
    gets selected is sensitive to numerical details of exactly how/when
    the perturbation lands, so this checks magnitude, not sign.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter24.ipynb")

    tail = ns.diff[-5:]
    assert (abs(tail) > 1.0).all(), f"expected a persistent ~1.5ms offset, got {tail}"
    assert (abs(tail) < 2.0).all(), f"offset magnitude too large: {tail}"
