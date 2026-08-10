from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_DIR = "27_Phase_Locking_with_Delays/THREE_DELAYED_PULSE_COUPLED_OSC"


def test_three_delayed_pulse_coupled_osc_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural/statistical properties (each
    # oscillator spikes roughly once per period, near-synchrony emerges
    # for the larger delay) rather than exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")

    assert set(py.i_spikes_1) == {0, 1, 2}
    assert set(py.i_spikes_2) == {0, 1, 2}

    counts_1 = np.bincount(py.i_spikes_1)
    counts_2 = np.bincount(py.i_spikes_2)
    assert np.max(counts_1) - np.min(counts_1) <= 2
    assert np.max(counts_2) - np.min(counts_2) <= 2

    # delta=0.55 run should be much closer to synchrony (spikes of the 3
    # oscillators clustered together) than the delta=0.45 run.
    def cluster_spread(t_spikes, i_spikes):
        # non-overlapping consecutive triplets only -- an overlapping
        # sliding window spuriously inflates the spread when spikes come
        # in exactly-synchronous triplets (ties get reordered arbitrarily
        # by argsort, so a window straddling two triplets looks "mixed").
        order = np.argsort(t_spikes)
        t_sorted = t_spikes[order]
        i_sorted = i_spikes[order]
        n_triplets = len(t_sorted) // 3
        spreads = []
        for k in range(n_triplets):
            chunk_i = i_sorted[3 * k:3 * k + 3]
            chunk_t = t_sorted[3 * k:3 * k + 3]
            if set(chunk_i) == {0, 1, 2}:
                spreads.append(chunk_t.max() - chunk_t.min())
        return np.median(spreads)

    spread_1 = cluster_spread(py.t_spikes_1, py.i_spikes_1)
    spread_2 = cluster_spread(py.t_spikes_2, py.i_spikes_2)
    assert spread_2 < spread_1
