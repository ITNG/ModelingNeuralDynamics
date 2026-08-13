from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_wb_with_depressing_s_structure():
    # deterministic (no RNG) -- verified visually against matlab's
    # figure.pdf (exact match).
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter39.ipynb")
    result = ns.simulate_wb_with_depressing_s()
    assert result.num_spikes > 3
    # depressing synapse: p recovers toward 1 between spikes but each
    # spike knocks it down sharply
    assert result.p.min() < 0.2
    assert result.p.max() == 1.0
