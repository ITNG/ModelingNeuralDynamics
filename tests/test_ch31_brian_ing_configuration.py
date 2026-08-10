from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_all_ten_ing_configurations_are_exact():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    expected = {
        "ING_1": (0.00, 0.5, 1.00, 0.00, 1.00, False, "uniform_random", 124875),
        "ING_2": (0.03, 0.5, 1.00, 0.00, 1.00, False, "phase", 63806),
        "ING_3": (0.00, 0.5, 0.85, 0.00, 1.00, False, "phase", 63806),
        "ING_4": (0.00, 0.5, 0.85, 0.00, 1.00, True, "phase", 63806),
        "ING_5": (0.05, 0.5, 0.50, 0.00, 1.00, False, "phase", 63806),
        "ING_6": (0.05, 0.5, 0.50, 0.10, 0.05, False, "phase", 63806),
        "ING_7": (0.00, 0.5, 1.00, 0.00, 1.00, False, "phase", 63806),
        "ING_8": (0.00, 0.5, 1.00, 0.00, 0.05, False, "phase", 63806),
        "ING_9": (0.05, 0.5, 1.00, 0.00, 0.05, False, "phase", 63806),
        "ING_10": (0.05, 0.5, 1.00, 0.04, 0.05, False, "phase", 63806),
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
            cfg["seed"],
        )
        for name, cfg in ns.ING_CONFIGS.items()
    }
    assert actual == expected


def test_connectivity_normalization_and_symmetry():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    built = ns.build_ing_inputs(ns.ING_CONFIGS["ING_6"], seed=63806, num_i=40)
    repeated = ns.build_ing_inputs(ns.ING_CONFIGS["ING_6"], seed=63806, num_i=40)
    assert built["g_ii_ms"].shape == (40, 40)
    assert built["g_gap_ms"].shape == (40, 40)
    assert np.allclose(built["g_gap_ms"], built["g_gap_ms"].T)
    assert np.allclose(np.diag(built["g_gap_ms"]), 0.0)
    assert np.array_equal(
        np.unique(built["g_ii_ms"][built["g_ii_ms"] != 0]),
        np.array([0.5 / (40 * 0.5)]),
    )
    assert np.array_equal(
        np.unique(built["g_gap_ms"][built["g_gap_ms"] != 0]),
        np.array([0.1 / (0.05 * 39)]),
    )
    for key in ("i_ext_ua", "g_ii_ms", "g_gap_ms"):
        assert np.array_equal(built[key], repeated[key])
    for key in ("v", "h", "n", "q", "s"):
        assert np.array_equal(
            built["initial_state"][key], repeated["initial_state"][key]
        )

    np.testing.assert_allclose(
        built["i_ext_ua"][:5],
        [
            1.5612064206593477,
            1.455691383444681,
            1.3904953729312304,
            1.3667372526218955,
            1.5796953639516627,
        ],
        rtol=0.0,
        atol=1e-14,
    )
    assert np.flatnonzero(built["g_ii_ms"][0]).tolist() == [
        0,
        3,
        7,
        9,
        10,
        12,
        14,
        16,
        18,
        19,
        21,
        22,
        23,
        24,
        27,
        31,
        34,
        35,
        36,
        37,
    ]
    assert np.argwhere(np.triu(built["g_gap_ms"], 1) != 0)[:10].tolist() == [
        [0, 3],
        [0, 6],
        [1, 29],
        [2, 6],
        [2, 7],
        [2, 17],
        [3, 10],
        [4, 13],
        [4, 22],
        [5, 7],
    ]
    np.testing.assert_allclose(
        built["initial_state"]["v"][:3] / ns.b2.mV,
        [-63.99251714811138, -60.465813304181275, -53.401452783580204],
        rtol=0.0,
        atol=1e-12,
    )


def test_fixed_indegree_has_equal_column_counts():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    built = ns.build_ing_inputs(ns.ING_CONFIGS["ING_4"], seed=63806, num_i=40)
    repeated = ns.build_ing_inputs(ns.ING_CONFIGS["ING_4"], seed=63806, num_i=40)
    counts = np.count_nonzero(built["g_ii_ms"], axis=0)
    assert np.array_equal(counts, np.full(40, round(0.85 * 40)))
    assert np.array_equal(built["g_ii_ms"], repeated["g_ii_ms"])
    assert np.flatnonzero(built["g_ii_ms"][:, 0] == 0).tolist() == [
        0,
        2,
        7,
        15,
        29,
        34,
    ]
    assert np.flatnonzero(built["g_ii_ms"][:, 1] == 0).tolist() == [
        2,
        13,
        23,
        24,
        26,
        31,
    ]
    assert np.flatnonzero(built["g_ii_ms"][:, 2] == 0).tolist() == [
        15,
        21,
        25,
        28,
        30,
        33,
    ]
    assert not np.array_equal(
        np.count_nonzero(built["g_ii_ms"], axis=1),
        np.full(40, round(0.85 * 40)),
    )
    assert np.array_equal(
        np.unique(built["g_ii_ms"][built["g_ii_ms"] != 0]),
        np.array([0.5 / (40 * 0.85)]),
    )


def test_zero_probability_normalization_is_explicit():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    config = dict(ns.ING_CONFIGS["ING_6"], p_ii=0.0)
    with pytest.raises(ValueError, match="g_hat_ii"):
        ns.build_ing_inputs(config, num_i=4)
    config = dict(ns.ING_CONFIGS["ING_6"], p_gap=0.0)
    with pytest.raises(ValueError, match="g_hat_gap"):
        ns.build_ing_inputs(config, num_i=4)

    config = dict(
        ns.ING_CONFIGS["ING_6"],
        g_hat_ii=0.0,
        p_ii=0.0,
        g_hat_gap=0.0,
        p_gap=0.0,
        init_mode="uniform_random",
    )
    built = ns.build_ing_inputs(config, seed=63806, num_i=4)
    assert np.array_equal(built["g_ii_ms"], np.zeros((4, 4)))
    assert np.array_equal(built["g_gap_ms"], np.zeros((4, 4)))


def test_wb_phase_initial_state_matches_second_to_third_crossing_oracle():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    state = ns.wb_phase_initial_state(
        np.full(3, 1.5), np.array([0.0, 0.5, 1.0])
    )
    assert set(state) == {"v", "h", "n", "q", "s"}
    assert all(value.shape == (3,) for value in state.values())
    np.testing.assert_allclose(
        state["v"] / ns.b2.mV,
        [-20.0, -60.207466427994085, -20.0],
        rtol=0.0,
        atol=1e-11,
    )
    np.testing.assert_allclose(
        state["h"],
        [0.022577870771095052, 0.7074439507194945, 0.022579384280192016],
        rtol=0.0,
        atol=1e-13,
    )
    np.testing.assert_allclose(
        state["n"],
        [0.7274002046819126, 0.10998835296871244, 0.72738877916241],
        rtol=0.0,
        atol=1e-13,
    )
    assert np.array_equal(state["q"], np.zeros(3))
    assert np.array_equal(state["s"], np.zeros(3))

    with pytest.raises(ValueError, match="phase"):
        ns.wb_phase_initial_state([1.5], [1.01])
    with pytest.raises(ValueError, match="equal-length"):
        ns.wb_phase_initial_state([1.5, 1.5], [0.5])
