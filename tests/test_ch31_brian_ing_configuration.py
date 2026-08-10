from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_all_ten_ing_configurations_are_exact():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    expected = {
        "ING_1": (0.00, 0.5, 1.00, 0.00, 1.00, False, "uniform_random"),
        "ING_2": (0.03, 0.5, 1.00, 0.00, 1.00, False, "phase"),
        "ING_3": (0.00, 0.5, 0.85, 0.00, 1.00, False, "phase"),
        "ING_4": (0.00, 0.5, 0.85, 0.00, 1.00, True, "phase"),
        "ING_5": (0.05, 0.5, 0.50, 0.00, 1.00, False, "phase"),
        "ING_6": (0.05, 0.5, 0.50, 0.10, 0.05, False, "phase"),
        "ING_7": (0.00, 0.5, 1.00, 0.00, 1.00, False, "phase"),
        "ING_8": (0.00, 0.5, 1.00, 0.00, 0.05, False, "phase"),
        "ING_9": (0.05, 0.5, 1.00, 0.00, 0.05, False, "phase"),
        "ING_10": (0.05, 0.5, 1.00, 0.04, 0.05, False, "phase"),
    }
    actual = {
        name: (
            cfg["sigma_i"],
            cfg["g_hat_ii"],
            cfg["p_ii"],
            cfg["g_hat_gap"],
            cfg["p_gap"],
            cfg["fixed_indegree"],
            cfg["init_mode"],
        )
        for name, cfg in ns.ING_CONFIGS.items()
    }
    assert actual == expected


def test_connectivity_normalization_and_symmetry():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    built = ns.build_ing_inputs(ns.ING_CONFIGS["ING_6"], seed=63806, num_i=40)
    assert built["g_ii_ms"].shape == (40, 40)
    assert built["g_gap_ms"].shape == (40, 40)
    assert np.allclose(built["g_gap_ms"], built["g_gap_ms"].T)
    assert np.allclose(np.diag(built["g_gap_ms"]), 0.0)


def test_fixed_indegree_has_equal_column_counts():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    built = ns.build_ing_inputs(ns.ING_CONFIGS["ING_4"], seed=63806, num_i=40)
    counts = np.count_nonzero(built["g_ii_ms"], axis=0)
    assert np.array_equal(counts, np.full(40, round(0.85 * 40)))
