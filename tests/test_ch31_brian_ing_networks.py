from pathlib import Path

import nbformat
import numpy as np
import pytest

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = ROOT / "python" / "31_ING_Rhythms"


class BuildCaptured(Exception):
    """Stop a simulator call immediately after capturing its input contract."""


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


def test_each_ing_name_selects_its_own_configuration_and_forwards_build_overrides(
    monkeypatch,
):
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    captured = {}

    def capture_build(config, seed=None, num_i=None):
        captured[current_name] = (dict(config), seed, num_i)
        raise BuildCaptured

    monkeypatch.setitem(
        ns.simulate_ing_network.__globals__, "build_ing_inputs", capture_build
    )
    for index in range(1, 11):
        current_name = f"ING_{index}"
        with pytest.raises(BuildCaptured):
            ns.simulate_ing_network(
                current_name,
                duration=0 * ns.b2.ms,
                num_i=10 + index,
                seed=900 + index,
            )

    assert list(captured) == [f"ING_{index}" for index in range(1, 11)]
    for index in range(1, 11):
        name = f"ING_{index}"
        config, seed, num_i = captured[name]
        assert config == ns.ING_CONFIGS[name]
        assert seed == 900 + index
        assert num_i == 10 + index


def test_synapse_probes_lock_matrix_orientation_weights_and_summed_targets(
    monkeypatch,
):
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    chemical_matrix = np.array(
        [
            [0.00, 0.11, 0.00],
            [0.22, 0.00, 0.23],
            [0.00, 0.00, 0.00],
        ]
    )
    gap_matrix = np.array(
        [
            [0.00, 0.31, 0.00],
            [0.32, 0.00, 0.33],
            [0.00, 0.34, 0.00],
        ]
    )

    def synthetic_inputs(config, seed=None, num_i=None):
        assert num_i == 3
        return {
            "i_ext_ua": np.array([1.4, 1.5, 1.6]),
            "g_ii_ms": chemical_matrix,
            "g_gap_ms": gap_matrix,
            "initial_state": {
                "v": np.array([-61.0, -62.0, -63.0]) * ns.b2.mV,
                "h": np.array([0.1, 0.2, 0.3]),
                "n": np.array([0.4, 0.5, 0.6]),
                "q": np.zeros(3),
                "s": np.zeros(3),
            },
        }

    probes = []

    class SynapsesProbe:
        def __init__(self, source, target, model):
            self.source = source
            self.target = target
            self.model = model
            self.i = None
            self.j = None
            self.g = None
            probes.append(self)

        def connect(self, i, j):
            self.i = np.asarray(i, dtype=int)
            self.j = np.asarray(j, dtype=int)

    monkeypatch.setitem(
        ns.simulate_ing_network.__globals__, "build_ing_inputs", synthetic_inputs
    )
    monkeypatch.setattr(ns.b2, "Synapses", SynapsesProbe)
    monkeypatch.setattr(ns.b2, "run", lambda duration: None)

    result = ns.simulate_ing_network(
        "ING_6", duration=0 * ns.b2.ms, num_i=3, record=False
    )

    assert len(probes) == 2
    chemical, gap = probes
    assert chemical.source is chemical.target
    assert chemical.model.splitlines() == [
        "g : siemens",
        "I_chem_post = g*s_pre*(-75*mV-v_post) : amp (summed)",
    ]
    np.testing.assert_array_equal(chemical.i, [0, 1, 1])
    np.testing.assert_array_equal(chemical.j, [1, 0, 2])
    np.testing.assert_allclose(chemical.g / ns.b2.msiemens, [0.11, 0.22, 0.23])

    assert gap.source is gap.target
    assert gap.model.splitlines() == [
        "g : siemens",
        "I_gap_post = g*(v_pre-v_post) : amp (summed)",
    ]
    np.testing.assert_array_equal(gap.i, [1, 0, 2, 1])
    np.testing.assert_array_equal(gap.j, [0, 1, 1, 2])
    np.testing.assert_allclose(gap.g / ns.b2.msiemens, [0.31, 0.32, 0.33, 0.34])

    assert set(result) == {"t_i_ms", "i_i", "duration_ms", "num_i", "config"}
    assert type(result["t_i_ms"]) is np.ndarray
    assert type(result["i_i"]) is np.ndarray
    assert np.issubdtype(result["t_i_ms"].dtype, np.floating)
    assert np.issubdtype(result["i_i"].dtype, np.integer)
    assert "state" not in result


def test_short_recorded_run_honors_duration_population_and_dt_overrides():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    ns.b2.prefs.codegen.target = "numpy"

    result = ns.simulate_ing_network(
        "ING_1",
        duration=0.2 * ns.b2.ms,
        dt=0.05 * ns.b2.ms,
        num_i=4,
        seed=17,
        record=True,
    )

    assert set(result) == {
        "t_i_ms",
        "i_i",
        "duration_ms",
        "num_i",
        "config",
        "state",
    }
    assert type(result["t_i_ms"]) is np.ndarray
    assert type(result["i_i"]) is np.ndarray
    assert np.issubdtype(result["t_i_ms"].dtype, np.floating)
    assert np.issubdtype(result["i_i"].dtype, np.integer)
    assert type(result["duration_ms"]) is float
    assert type(result["num_i"]) is int
    assert type(result["config"]) is dict
    assert result["duration_ms"] == pytest.approx(0.2)
    assert result["num_i"] == 4
    assert result["config"] == ns.ING_CONFIGS["ING_1"]

    state = result["state"]
    assert set(state) == {"t_ms", "v_mv", "q", "s"}
    assert all(type(values) is np.ndarray for values in state.values())
    np.testing.assert_allclose(state["t_ms"], [0.0, 0.05, 0.10, 0.15])
    assert state["v_mv"].shape == (4, 4)
    assert state["q"].shape == (4, 4)
    assert state["s"].shape == (4, 4)


def test_unknown_ing_name_is_rejected_before_start_scope(monkeypatch):
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    start_scope_calls = []

    def capture_start_scope():
        start_scope_calls.append(True)

    monkeypatch.setattr(ns.b2, "start_scope", capture_start_scope)
    with pytest.raises(ValueError, match="ING_99"):
        ns.simulate_ing_network("ING_99")
    assert start_scope_calls == []


def test_plot_ing_raster_draws_spikes_and_full_population_bounds():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    result = {
        "t_i_ms": np.array([1.0, 3.0, 7.0]),
        "i_i": np.array([2, 1, 3]),
        "duration_ms": 20.0,
        "num_i": 4,
    }

    ax = ns.plot_ing_raster(result)

    assert len(ax.lines) == 1
    np.testing.assert_array_equal(ax.lines[0].get_xdata(), [1.0, 3.0, 7.0])
    np.testing.assert_array_equal(ax.lines[0].get_ydata(), [2, 1, 3])
    assert ax.lines[0].get_marker() == "."
    assert ax.lines[0].get_color() == "k"
    assert ax.lines[0].get_markersize() == 2.0
    assert tuple(ax.get_xlim()) == (0.0, 20.0)
    assert tuple(ax.get_ylim()) == (0.0, 5.0)
    assert ax.get_xlabel() == "t [ms]"
    ns.plt.close(ax.figure)


def test_plot_ing_raster_reuses_axis_and_applies_closeup_bounds():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    result = {
        "t_i_ms": np.array([410.0, 450.0, 490.0]),
        "i_i": np.array([1, 10, 20]),
        "duration_ms": 500.0,
        "num_i": 100,
    }
    _, supplied_ax = ns.plt.subplots()

    returned_ax = ns.plot_ing_raster(result, ax=supplied_ax, closeup=True)

    assert returned_ax is supplied_ax
    assert tuple(supplied_ax.get_xlim()) == (400.0, 500.0)
    assert tuple(supplied_ax.get_ylim()) == (0.0, 21.0)
    ns.plt.close(supplied_ax.figure)


def test_guarded_raster_cells_run_all_variants_with_required_closeups():
    notebook = nbformat.read(ROOT / "brian" / "chapter31.ipynb", as_version=4)
    raster_cells = [
        cell
        for cell in notebook.cells
        if cell.cell_type == "code" and cell.id.endswith("-raster")
    ]
    assert len(raster_cells) == 10
    source = "\n".join("".join(cell.source) for cell in raster_cells)
    calls = []

    def capture_simulation(name):
        calls.append(("simulate", name))
        return {"name": name}

    def capture_plot(result, ax=None, closeup=False):
        calls.append(("plot", result["name"], closeup))

    namespace = {
        "__name__": "__main__",
        "simulate_ing_network": capture_simulation,
        "plot_ing_raster": capture_plot,
    }
    exec(compile(source, str(ROOT / "brian" / "chapter31.ipynb"), "exec"), namespace)

    expected_names = [f"ING_{index}" for index in range(1, 11)]
    assert calls[::2] == [("simulate", name) for name in expected_names]
    assert calls[1::2] == [
        ("plot", name, name in {"ING_8", "ING_9", "ING_10"})
        for name in expected_names
    ]

    calls.clear()
    namespace["__name__"] = "chapter31"
    exec(compile(source, str(ROOT / "brian" / "chapter31.ipynb"), "exec"), namespace)
    assert calls == []
