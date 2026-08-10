from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = ROOT / "python" / "31_ING_Rhythms"


@pytest.mark.parametrize("name", [f"ING_{index}" for index in range(1, 11)])
def test_brian_ing_network_matches_python_population_activity(name):
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    result = ns.simulate_ing_network(name)
    assert len(result["t_i_ms"]) > 500
    counts = np.bincount(result["i_i"].astype(int), minlength=result["num_i"])
    assert np.count_nonzero(counts) > 0.9 * result["num_i"]

    if name != "ING_1":
        py = load_python_port(PYTHON_BASE / name / "main.py")
        brian_rate = len(result["t_i_ms"]) / (
            result["num_i"] * result["duration_ms"] / 1000.0
        )
        python_rate = len(py.t_i_spikes) / (py.num_i * py.t_final / 1000.0)
        assert np.isclose(brian_rate, python_rate, rtol=0.20)


def test_ing_network_rejects_unknown_name():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    with pytest.raises(ValueError, match="ING_99"):
        ns.simulate_ing_network("ING_99", duration=20 * ns.b2.ms, num_i=10)
