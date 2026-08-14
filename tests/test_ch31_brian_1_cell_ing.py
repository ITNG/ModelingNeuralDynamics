from pathlib import Path

import numpy as np

from matlab_ref import (
    load_notebook_as_module,
    load_notebook_definitions_as_module,
    spike_times,
    trace_rmse,
)

ROOT = Path(__file__).resolve().parents[1]


def test_brian_single_cell_ing_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    py = load_notebook_definitions_as_module(ROOT / "python" / "chapter31.ipynb")
    t, v = py.simulate_1_cell_ing()
    expected_spikes = spike_times(t, v, threshold=-20.0)

    result = ns.simulate_single_cell_ing()
    actual_spikes = np.asarray(result["spike_ms"])

    assert trace_rmse(t, v, result["t_ms"], result["v_mv"]) < 3.0
    assert len(actual_spikes) == len(expected_spikes)
    assert np.isclose(result["period_ms"], np.diff(expected_spikes)[-1], atol=0.5)
