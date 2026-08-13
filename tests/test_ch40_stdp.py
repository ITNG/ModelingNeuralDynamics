from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_three_cell_ping_5_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy (though this particular script draws no random numbers, so
    # results are actually deterministic); verified visually against
    # matlab's figure.pdf (close match: g_12 decays toward 0, g_21
    # grows toward its bound B, and the E1-E2 lag shrinks over time).
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")
    py = ns.simulate_three_cell_ping_5()
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
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")
    py = ns.simulate_ping_with_stdp()
    assert len(py.t_e_spikes) > 500
    assert len(py.t_i_spikes) > 100
    assert py.vec_g.min() >= 0.0
    assert py.vec_g.max() <= py.B.max() * 1.01
    # STDP should have produced a spread of synaptic strengths, not
    # left every synapse exactly at its initial value
    assert py.vec_g.std() > 0
