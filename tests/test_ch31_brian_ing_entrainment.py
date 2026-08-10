from pathlib import Path

import nbformat
import numpy as np
import pytest

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = ROOT / "python" / "31_ING_Rhythms"


def _active_fraction(indices, count):
    return np.count_nonzero(
        np.bincount(np.asarray(indices, dtype=int), minlength=count)
    ) / count


@pytest.fixture(scope="module")
def ns():
    return load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")


def _zero_state(ns, count):
    return {
        "v": np.full(count, -70.0) * ns.b2.mV,
        "h": np.zeros(count),
        "n": np.zeros(count),
        "q": np.zeros(count),
        "s": np.zeros(count),
    }


def test_entrainment_configurations_are_exact(ns):
    assert ns.ENTRAINMENT_CONFIGS == {
        "ING_ENTRAINING_E_CELLS": {
            "num_e": 400,
            "num_i": 100,
            "sigma_e": 0.10,
            "sigma_i": 0.05,
            "mean_e": 1.5,
            "mean_i": 1.5,
            "g_hat_ee": 0.0,
            "g_hat_ei": 0.0,
            "g_hat_ie": 0.50,
            "g_hat_ii": 0.50,
            "p_ee": 1.0,
            "p_ei": 0.5,
            "p_ie": 0.5,
            "p_ii": 0.5,
            "g_hat_gap": 0.1,
            "p_gap": 0.05,
        },
        "ING_ENTRAINING_E_CELLS_2": {
            "num_e": 400,
            "num_i": 100,
            "sigma_e": 0.0,
            "sigma_i": 0.0,
            "mean_e": 1.9,
            "mean_i": 1.5,
            "g_hat_ee": 0.0,
            "g_hat_ei": 0.0,
            "g_hat_ie": 0.50,
            "g_hat_ii": 0.50,
            "p_ee": 1.0,
            "p_ei": 1.0,
            "p_ie": 1.0,
            "p_ii": 1.0,
            "g_hat_gap": 0.1,
            "p_gap": 0.05,
        },
    }


def test_task6_wb_phase_origin_matches_reference_oracle(ns):
    state = ns.wb_entrainment_phase_initial_state(
        np.array([1.45, 1.50, 1.55]), np.array([0.2, 0.5, 0.8])
    )
    actual = np.column_stack(
        (state["v"] / ns.b2.mV, state["h"], state["n"])
    )
    expected = np.array(
        [
            [-65.28114058195158, 0.6097889310143990, 0.1632596427880604],
            [-60.20724162773239, 0.7074370853511656, 0.10999015164876737],
            [-54.05019224159273, 0.5588706808553839, 0.15639186573826044],
        ]
    )
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=3e-11)
    assert np.array_equal(state["q"], np.zeros(3))
    assert np.array_equal(state["s"], np.zeros(3))


def test_task6_rtm_phase_origin_matches_reference_oracle(ns):
    state = ns.rtm_phase_initial_state(
        np.array([1.5, 1.9, 2.1]), np.array([0.2, 0.5, 0.8])
    )
    actual = np.column_stack(
        (state["v"] / ns.b2.mV, state["h"], state["n"])
    )
    expected = np.array(
        [
            [-86.43789772277115, 0.9926543238355179, 0.01082429435913847],
            [-72.93397154144637, 0.9990074542870481, 0.008779755975138186],
            [-63.55474270688976, 0.9962048349605402, 0.04302880085010221],
        ]
    )
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=3e-11)
    assert np.array_equal(state["q"], np.zeros(3))
    assert np.array_equal(state["s"], np.zeros(3))


def test_build_entrainment_inputs_locks_rng_order_and_normalization(ns, monkeypatch):
    calls = []

    def capture_rtm(i_ext, phase):
        calls.append(("rtm", np.asarray(i_ext).copy(), np.asarray(phase).copy()))
        return _zero_state(ns, len(i_ext))

    def capture_wb(i_ext, phase):
        calls.append(("wb", np.asarray(i_ext).copy(), np.asarray(phase).copy()))
        return _zero_state(ns, len(i_ext))

    def reject_legacy_wb(*args, **kwargs):
        raise AssertionError("Task 6 must not use the second-to-third WB initializer")

    globals_ = ns.build_entrainment_inputs.__globals__
    monkeypatch.setitem(globals_, "rtm_phase_initial_state", capture_rtm)
    monkeypatch.setitem(globals_, "wb_entrainment_phase_initial_state", capture_wb)
    monkeypatch.setitem(globals_, "wb_phase_initial_state", reject_legacy_wb)

    built = ns.build_entrainment_inputs(
        ns.ENTRAINMENT_CONFIGS["ING_ENTRAINING_E_CELLS"],
        seed=17,
        num_e=3,
        num_i=3,
    )

    np.testing.assert_allclose(
        built["i_ext_e_ua"],
        [1.6651893680258771, 1.5507646914969269, 1.4190042727119745],
        rtol=0.0,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        built["i_ext_i_ua"],
        [1.4054818607360675, 1.3579034047620557, 1.5013978718237464],
        rtol=0.0,
        atol=1e-15,
    )
    np.testing.assert_array_equal(built["g_ee_ms"], np.zeros((3, 3)))
    np.testing.assert_array_equal(built["g_ei_ms"], np.zeros((3, 3)))
    np.testing.assert_array_equal(
        built["g_ie_ms"],
        np.array(
            [
                [0.0, 1 / 3, 1 / 3],
                [1 / 3, 0.0, 1 / 3],
                [0.0, 0.0, 1 / 3],
            ]
        ),
    )
    np.testing.assert_array_equal(
        built["g_ii_ms"],
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1 / 3, 0.0, 1 / 3],
                [0.0, 0.0, 1 / 3],
            ]
        ),
    )
    np.testing.assert_array_equal(
        built["g_gap_ms"],
        np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
    )
    assert [kind for kind, _, _ in calls] == ["rtm", "wb"]
    np.testing.assert_allclose(
        calls[0][1], built["i_ext_e_ua"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        calls[1][1], built["i_ext_i_ua"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        calls[0][2],
        [0.1713635614270761, 0.3430499725407290, 0.6175147346249099],
        rtol=0.0,
        atol=1e-15,
    )
    np.testing.assert_allclose(
        calls[1][2],
        [0.5508373261068352, 0.7160952791161002, 0.4894537339227958],
        rtol=0.0,
        atol=1e-15,
    )


def test_all_four_chemical_normalizations_use_presynaptic_population(ns, monkeypatch):
    monkeypatch.setitem(
        ns.build_entrainment_inputs.__globals__,
        "rtm_phase_initial_state",
        lambda i_ext, phase: _zero_state(ns, len(i_ext)),
    )
    monkeypatch.setitem(
        ns.build_entrainment_inputs.__globals__,
        "wb_entrainment_phase_initial_state",
        lambda i_ext, phase: _zero_state(ns, len(i_ext)),
    )
    config = {
        **ns.ENTRAINMENT_CONFIGS["ING_ENTRAINING_E_CELLS_2"],
        "g_hat_ee": 0.2,
        "g_hat_ei": 0.3,
        "g_hat_ie": 0.4,
        "g_hat_ii": 0.5,
        "g_hat_gap": 0.0,
        "p_gap": 1.0,
    }
    built = ns.build_entrainment_inputs(config, seed=5, num_e=2, num_i=3)
    np.testing.assert_allclose(built["g_ee_ms"], np.full((2, 2), 0.2 / 2))
    np.testing.assert_allclose(built["g_ei_ms"], np.full((2, 3), 0.3 / 2))
    np.testing.assert_allclose(built["g_ie_ms"], np.full((3, 2), 0.4 / 3))
    np.testing.assert_allclose(built["g_ii_ms"], np.full((3, 3), 0.5 / 3))


class SynapsesProbe:
    instances = []

    def __init__(self, source, target, model):
        self.source = source
        self.target = target
        self.model = model
        self.i = None
        self.j = None
        self.condition = None
        self.g = None
        self.E_rev = None
        self.instances.append(self)

    def connect(self, i=None, j=None, condition=None):
        self.i = None if i is None else np.asarray(i, dtype=int)
        self.j = None if j is None else np.asarray(j, dtype=int)
        self.condition = condition


def test_simulator_locks_projection_contract_forwarding_and_result_types(ns, monkeypatch):
    matrices = {
        "g_ee_ms": np.array([[0.0, 0.11], [0.12, 0.0]]),
        "g_ei_ms": np.array([[0.21, 0.0, 0.23], [0.0, 0.24, 0.0]]),
        "g_ie_ms": np.array([[0.31, 0.0], [0.0, 0.32], [0.33, 0.34]]),
        "g_ii_ms": np.array(
            [[0.0, 0.41, 0.0], [0.42, 0.0, 0.43], [0.0, 0.44, 0.0]]
        ),
        "g_gap_ms": np.array(
            [[0.0, 0.51, 0.0], [0.52, 0.0, 0.53], [0.0, 0.54, 0.0]]
        ),
    }
    build_calls = []
    run_calls = []

    def synthetic_inputs(config, seed, num_e, num_i):
        build_calls.append((dict(config), seed, num_e, num_i))
        return {
            "i_ext_e_ua": np.array([1.9, 2.0]),
            "i_ext_i_ua": np.array([1.4, 1.5, 1.6]),
            **matrices,
            "initial_e": _zero_state(ns, 2),
            "initial_i": _zero_state(ns, 3),
        }

    SynapsesProbe.instances = []
    monkeypatch.setitem(
        ns.simulate_ing_entrainment.__globals__,
        "build_entrainment_inputs",
        synthetic_inputs,
    )
    monkeypatch.setattr(ns.b2, "Synapses", SynapsesProbe)
    monkeypatch.setattr(ns.b2, "run", lambda duration: run_calls.append(duration))

    result = ns.simulate_ing_entrainment(
        "ING_ENTRAINING_E_CELLS_2",
        e_drive=2.0,
        duration=0.3 * ns.b2.ms,
        dt=0.05 * ns.b2.ms,
        num_e=2,
        num_i=3,
        seed=41,
    )

    assert len(build_calls) == 1
    config, seed, num_e, num_i = build_calls[0]
    assert config == {
        **ns.ENTRAINMENT_CONFIGS["ING_ENTRAINING_E_CELLS_2"],
        "mean_e": 2.0,
        "num_e": 2,
        "num_i": 3,
    }
    assert (seed, num_e, num_i) == (41, 2, 3)
    assert len(run_calls) == 1
    assert float(run_calls[0] / ns.b2.ms) == pytest.approx(0.3)
    assert float(ns.b2.defaultclock.dt / ns.b2.ms) == pytest.approx(0.05)

    ee, ei, ie, ii, gap = SynapsesProbe.instances
    excitatory = ee.source
    inhibitory = ei.target
    assert (ee.source, ee.target) == (excitatory, excitatory)
    assert (ei.source, ei.target) == (excitatory, inhibitory)
    assert (ie.source, ie.target) == (inhibitory, excitatory)
    assert (ii.source, ii.target) == (inhibitory, inhibitory)
    assert (gap.source, gap.target) == (inhibitory, inhibitory)

    for probe, target_variable, reversal in (
        (ee, "I_from_e", 0.0),
        (ei, "I_from_e", 0.0),
        (ie, "I_from_i", -75.0),
        (ii, "I_from_i", -75.0),
    ):
        assert probe.model.splitlines() == [
            "g : siemens",
            "E_rev : volt (constant)",
            f"{target_variable}_post = g*s_pre*(E_rev-v_post) : amp (summed)",
        ]
        assert float(probe.E_rev / ns.b2.mV) == pytest.approx(reversal)

    expected_connections = (
        ([0, 1], [1, 0], [0.11, 0.12]),
        ([0, 0, 1], [0, 2, 1], [0.21, 0.23, 0.24]),
        ([0, 1, 2, 2], [0, 1, 0, 1], [0.31, 0.32, 0.33, 0.34]),
        ([0, 1, 1, 2], [1, 0, 2, 1], [0.41, 0.42, 0.43, 0.44]),
    )
    for probe, (pre, post, weights) in zip((ee, ei, ie, ii), expected_connections):
        np.testing.assert_array_equal(probe.i, pre)
        np.testing.assert_array_equal(probe.j, post)
        np.testing.assert_allclose(probe.g / ns.b2.msiemens, weights)

    assert gap.model.splitlines() == [
        "g : siemens",
        "I_gap_post = g*(v_pre-v_post) : amp (summed)",
    ]
    np.testing.assert_array_equal(gap.i, [1, 0, 2, 1])
    np.testing.assert_array_equal(gap.j, [0, 1, 1, 2])
    np.testing.assert_allclose(
        gap.g / ns.b2.msiemens, [0.51, 0.52, 0.53, 0.54]
    )

    assert set(result) == {
        "t_e_ms",
        "i_e",
        "t_i_ms",
        "i_i",
        "duration_ms",
        "num_e",
        "num_i",
        "config",
    }
    for key in ("t_e_ms", "i_e", "t_i_ms", "i_i"):
        assert type(result[key]) is np.ndarray
    assert np.issubdtype(result["t_e_ms"].dtype, np.floating)
    assert np.issubdtype(result["t_i_ms"].dtype, np.floating)
    assert np.issubdtype(result["i_e"].dtype, np.integer)
    assert np.issubdtype(result["i_i"].dtype, np.integer)
    assert type(result["duration_ms"]) is float
    assert type(result["num_e"]) is int
    assert type(result["num_i"]) is int
    assert type(result["config"]) is dict
    assert result["duration_ms"] == pytest.approx(0.3)
    assert result["num_e"] == 2
    assert result["num_i"] == 3


@pytest.mark.parametrize(
    ("variant", "kwargs", "match"),
    [
        ("UNKNOWN", {}, "UNKNOWN"),
        ("ING_ENTRAINING_E_CELLS", {"num_e": 0}, "num_e"),
        ("ING_ENTRAINING_E_CELLS", {"num_i": 0}, "num_i"),
        ("ING_ENTRAINING_E_CELLS", {"e_drive": 2.0}, "only selectable"),
        ("ING_ENTRAINING_E_CELLS_2", {"e_drive": 2.2}, "one of"),
    ],
)
def test_invalid_entrainment_requests_fail_before_start_scope(
    ns, monkeypatch, variant, kwargs, match
):
    start_scope_calls = []
    monkeypatch.setattr(ns.b2, "start_scope", lambda: start_scope_calls.append(True))
    call_kwargs = dict(kwargs)
    with pytest.raises(ValueError, match=match):
        ns.simulate_ing_entrainment(
            variant, duration=0 * ns.b2.ms, num_e=call_kwargs.pop("num_e", 2),
            num_i=call_kwargs.pop("num_i", 2), **call_kwargs
        )
    assert start_scope_calls == []


def _plot_result(ns, drive=1.9):
    return {
        "t_e_ms": np.array([2.0, 8.0]),
        "i_e": np.array([1, 3]),
        "t_i_ms": np.array([1.0, 7.0]),
        "i_i": np.array([2, 4]),
        "duration_ms": 20.0,
        "num_e": 5,
        "num_i": 6,
        "config": {"mean_e": drive},
    }


def test_entrainment_plots_lock_colors_offsets_limits_and_labels(ns):
    result = _plot_result(ns)
    ax = ns.plot_ing_entrainment(result)
    assert len(ax.lines) == 3
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), [1.0, 7.0])
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), [2, 4])
    assert ax.lines[0].get_color() == "b"
    np.testing.assert_array_equal(ax.lines[1].get_xdata(), [2.0, 8.0])
    np.testing.assert_array_equal(ax.lines[1].get_ydata(), [7, 9])
    assert ax.lines[1].get_color() == "r"
    np.testing.assert_array_equal(ax.lines[2].get_ydata(), [6.5, 6.5])
    assert ax.lines[2].get_linestyle() == "--"
    assert ax.lines[2].get_color() == "k"
    assert tuple(ax.get_xlim()) == (0.0, 20.0)
    assert tuple(ax.get_ylim()) == (0.0, 12.0)
    assert ax.get_xlabel() == "t [ms]"
    ns.plt.close(ax.figure)

    _, axes = ns.plt.subplots(3, 1)
    returned = ns.plot_ing_entrainment_sweep(
        [_plot_result(ns, drive) for drive in (1.9, 2.0, 2.1)], axes=axes
    )
    assert returned.tolist() == axes.tolist()
    assert [ax.get_title() for ax in axes] == [
        r"$\overline{I}_E=1.9$",
        r"$\overline{I}_E=2$",
        r"$\overline{I}_E=2.1$",
    ]
    assert [ax.get_xlabel() for ax in axes] == ["", "", "t [ms]"]
    ns.plt.close(axes[0].figure)


def test_entrainment_plot_cells_are_guarded_and_run_expected_variants():
    notebook = nbformat.read(ROOT / "brian" / "chapter31.ipynb", as_version=4)
    cells = [
        cell
        for cell in notebook.cells
        if cell.id in {"ing-entraining-e-cells-plot", "ing-entraining-e-cells-2-plot"}
    ]
    assert [cell.id for cell in cells] == [
        "ing-entraining-e-cells-plot",
        "ing-entraining-e-cells-2-plot",
    ]
    source = "\n".join("".join(cell.source) for cell in cells)
    calls = []

    def capture_simulation(variant, e_drive=None):
        calls.append(("simulate", variant, e_drive))
        return {"variant": variant, "e_drive": e_drive}

    def capture_plot(result):
        calls.append(("plot", result["variant"], result["e_drive"]))

    def capture_sweep(results):
        calls.append(("sweep", tuple(result["e_drive"] for result in results)))

    namespace = {
        "__name__": "__main__",
        "simulate_ing_entrainment": capture_simulation,
        "plot_ing_entrainment": capture_plot,
        "plot_ing_entrainment_sweep": capture_sweep,
    }
    exec(compile(source, str(ROOT / "brian" / "chapter31.ipynb"), "exec"), namespace)
    assert calls == [
        ("simulate", "ING_ENTRAINING_E_CELLS", None),
        ("plot", "ING_ENTRAINING_E_CELLS", None),
        ("simulate", "ING_ENTRAINING_E_CELLS_2", 1.9),
        ("simulate", "ING_ENTRAINING_E_CELLS_2", 2.0),
        ("simulate", "ING_ENTRAINING_E_CELLS_2", 2.1),
        ("sweep", (1.9, 2.0, 2.1)),
    ]

    calls.clear()
    namespace["__name__"] = "chapter31"
    exec(compile(source, str(ROOT / "brian" / "chapter31.ipynb"), "exec"), namespace)
    assert calls == []


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
