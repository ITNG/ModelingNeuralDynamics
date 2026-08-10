from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = ROOT / "python" / "31_ING_Rhythms"


def _active_fraction(indices, count):
    return np.count_nonzero(
        np.bincount(np.asarray(indices, dtype=int), minlength=count)
    ) / count


def test_brian_ing_entraining_e_cells_matches_python_structure():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    py = load_python_port(PYTHON_BASE / "ING_ENTRAINING_E_CELLS" / "main.py")
    result = ns.simulate_ing_entrainment("ING_ENTRAINING_E_CELLS")
    assert _active_fraction(result["i_e"], result["num_e"]) > 0.9
    assert _active_fraction(result["i_i"], result["num_i"]) > 0.9
    assert np.isclose(len(result["t_e_ms"]), len(py.t_e_spikes), rtol=0.25)
    assert np.isclose(len(result["t_i_ms"]), len(py.t_i_spikes), rtol=0.25)


def test_brian_ing_entraining_e_cells_drive_sweep_has_three_results():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    results = [
        ns.simulate_ing_entrainment("ING_ENTRAINING_E_CELLS_2", e_drive=drive)
        for drive in (1.9, 2.0, 2.1)
    ]
    assert len(results) == 3
    assert all(len(result["t_e_ms"]) > 500 for result in results)
    assert all(len(result["t_i_ms"]) > 200 for result in results)
