from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_ping_6_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing in
    # all three connectivity panels) rather than exact spike times.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter30.ipynb")
    panels = ns.run_connectivity_panels(t_final_run=200.0)
    assert len(panels) == 3
    for t_e, i_e, t_i, i_i in panels:
        assert len(t_e) > 1000
        assert len(t_i) > 200
        assert len(t_e) == len(i_e)
        assert len(t_i) == len(i_i)


def test_ping_6_default_main_preserves_reference_plot(monkeypatch):
    # PING_6's main() has no `if __name__ == "__main__":` guard -- unlike
    # the legacy per-example main.py, the notebook's guard against
    # accidentally running heavy plotting code on import is provided by
    # the harness itself: load_notebook_definitions_as_module only ever
    # executes top-level imports/def statements from the notebook cells,
    # never top-level calls, so `main()` (like every other example's call
    # cell) only runs when a caller explicitly invokes it.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter30.ipynb")
    num_e, num_i, t_final = 200, 50, 2000.0

    calls = []
    panel = (
        np.array([10.0]),
        np.array([1]),
        np.array([20.0]),
        np.array([2]),
    )

    def fake_run_connectivity_panels():
        calls.append(())
        return [panel, panel, panel]

    monkeypatch.setitem(
        ns.main.__globals__, "run_connectivity_panels", fake_run_connectivity_panels
    )
    saved = []
    monkeypatch.setattr(ns.plt, "savefig", saved.append)

    ns.main()
    fig = ns.plt.gcf()
    try:
        assert calls == [()]
        assert saved == ["fig.png"]
        assert len(fig.axes) == 3
        for ax in fig.axes:
            np.testing.assert_array_equal(ax.get_yticks(), [num_i, num_e + num_i])
            assert any(
                np.array_equal(line.get_ydata(), [num_i + 0.5, num_i + 0.5])
                for line in ax.lines
            )
            np.testing.assert_allclose(
                ax.axis(),
                [t_final - 200, t_final, 0, num_e + num_i + 1],
            )
        assert fig.axes[-1].get_xlabel() == '$t$ [ms]'
    finally:
        ns.plt.close(fig)


def test_ping_7_structure():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter30.ipynb")
    t_e, i_e, t_i, i_i, lfp, num_e, num_i, t_final = ns.simulate_ping_7()
    assert len(t_e) > 200
    assert len(t_i) > 30
    assert lfp.min() > -100 and lfp.max() < 50


def test_ping_8_structure():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter30.ipynb")
    t_e, i_e, t_i, i_i, lfp, num_e, num_i, t_final = ns.simulate_ping_8()
    assert len(t_e) > 200
    assert len(t_i) > 30
    assert lfp.min() > -100 and lfp.max() < 50


def test_ping_9_structure():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter30.ipynb")
    t_e, i_e, t_i, i_i, lfp, num_e, num_i, t_final = ns.simulate_ping_9()
    assert len(t_e) > 200
    assert len(t_i) > 30
    assert lfp.min() > -100 and lfp.max() < 50
