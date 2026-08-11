from pathlib import Path

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_wilson_cowan_phase_plane_smoke():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter22.ipynb")
    ns.plot_wilson_cowan_phase_plane()
