from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_wilson_cowan_rastergram_smoke():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter22.ipynb")
    t, spikes_e, spikes_i = ns.simulate_wilson_cowan_rastergram(t_final=50.0)

    assert len(spikes_e) == 80
    assert len(spikes_i) == 20
