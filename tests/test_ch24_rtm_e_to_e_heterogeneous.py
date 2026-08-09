from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "24_Synchronization_by_Fast_Recurrent_Excitation/RTM_E_TO_E_HETEROGENEOUS"


def test_rtm_e_to_e_heterogeneous_is_sane():
    # matlab seeds its own RNG (rng(63806)) for the random i_ext/g_syn,
    # which numpy cannot reproduce bit-for-bit -- verified visually
    # against figure.pdf (same qualitative near-synchronous wave pattern,
    # same rough period) rather than by an exact numeric match here.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    assert py.i_ext.shape == (py.N,)
    assert np.all((0.25 <= py.i_ext) & (py.i_ext <= 0.35))

    assert py.g_syn.shape == (py.N, py.N)
    assert np.all(np.diag(py.g_syn) == 0)
    off_diag = py.g_syn[~np.eye(py.N, dtype=bool)]
    assert np.all((0.00625 <= off_diag) & (off_diag <= 0.00875))

    # every neuron should have spiked multiple times over the 200ms run
    counts = np.array([np.sum(py.i_spikes == i) for i in range(1, py.N + 1)])
    assert np.all(counts >= 3)
