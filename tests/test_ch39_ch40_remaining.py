from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
CH39_BASE = "39_Short-Term_Depression_and_Facilitation"
CH40_BASE = "40_Spike_Timing-Dependent_Plasticity(STDP)"


def test_wb_with_depressing_s_structure():
    # deterministic (no RNG) -- verified visually against matlab's
    # figure.pdf (exact match).
    py = load_python_port(ROOT / "python" / CH39_BASE / "WB_WITH_DEPRESSING_S" / "main.py")
    assert py.num_spikes > 3
    # depressing synapse: p recovers toward 1 between spikes but each
    # spike knocks it down sharply
    assert py.p.min() < 0.2
    assert py.p.max() == 1.0


def test_three_cell_ping_5_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy (though this particular script draws no random numbers, so
    # results are actually deterministic); verified visually against
    # matlab's figure.pdf (close match: g_12 decays toward 0, g_21
    # grows toward its bound B, and the E1-E2 lag shrinks over time).
    py = load_python_port(ROOT / "python" / CH40_BASE / "THREE_CELL_PING_5" / "main.py")
    assert len(py.t_e_spikes) > 20
    assert len(py.t_i_spikes) > 10
    assert py.g_12[-1] < py.g_12[0]
    assert py.g_21[-1] > py.g_21[0]
    assert py.g_21[-1] > 0.9 * py.B[1, 0]
    assert len(py.lags) > 5
    assert py.lags[-1] < py.lags[0]


def test_ping_with_stdp_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties rather than exact spikes.
    py = load_python_port(ROOT / "python" / CH40_BASE / "PING_WITH_STDP" / "main.py")
    assert len(py.t_e_spikes) > 500
    assert len(py.t_i_spikes) > 100
    assert py.vec_g.min() >= 0.0
    assert py.vec_g.max() <= py.B.max() * 1.01
    # STDP should have produced a spread of synaptic strengths, not
    # left every synapse exactly at its initial value
    assert py.vec_g.std() > 0
