from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]
NUM_E = 400
NUM_I = 100
T_FINAL = 500.0
I_EXT_E_VEC = (1.9, 2.0, 2.1)


def _load_ns():
    return load_notebook_definitions_as_module(ROOT / "python" / "chapter31.ipynb")


def test_ing_entraining_e_cells_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial E and I
    # firing, most cells active) rather than exact spike times.
    ns = _load_ns()
    t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, num_e, num_i = ns.simulate_ing_entraining_e_cells()
    assert len(t_e_spikes) > 1000
    assert len(t_i_spikes) > 200
    counts_e = np.bincount(i_e_spikes.astype(int), minlength=num_e)
    counts_i = np.bincount(i_i_spikes.astype(int), minlength=num_i)
    assert np.count_nonzero(counts_e) > 0.9 * num_e
    assert np.count_nonzero(counts_i) > 0.9 * num_i


def test_ing_entraining_e_cells_2_structure():
    ns = _load_ns()
    assert not hasattr(ns, "results")
    results = ns.run_drive_panels(t_final_run=100.0)
    assert len(results) == 3
    for t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes in results:
        assert len(t_e_spikes) > 500
        assert len(t_i_spikes) > 200
        assert len(t_e_spikes) == len(i_e_spikes)
        assert len(t_i_spikes) == len(i_i_spikes)


def test_ing_entraining_e_cells_2_default_main_preserves_reference_plot(monkeypatch):
    # Notebook-based port of the original main.py's "run_drive_panels then
    # render a fixed 3-panel figure" behavior. The original script also
    # guarded this behind `if __name__ == "__main__": main()`; that
    # guard's job (don't run the heavy simulation just from importing the
    # module) is now covered by load_notebook_definitions_as_module
    # itself, which only executes function/class/import statements and
    # never calls main() on its own.
    ns = _load_ns()
    calls = []
    panel = (
        np.array([10.0]),
        np.array([1]),
        np.array([20.0]),
        np.array([2]),
    )

    def fake_run_drive_panels():
        calls.append(())
        return [panel, panel, panel]

    monkeypatch.setitem(ns.main.__globals__, "run_drive_panels", fake_run_drive_panels)
    saved = []
    monkeypatch.setattr(ns.plt, "savefig", saved.append)

    ns.main()
    fig = ns.plt.gcf()
    try:
        assert calls == [()]
        assert saved == ["fig.png"]
        assert len(fig.axes) == 3
        for ax in fig.axes:
            np.testing.assert_array_equal(ax.get_yticks(), [NUM_I, NUM_E + NUM_I])
            assert any(
                np.array_equal(line.get_ydata(), [NUM_I + 0.5, NUM_I + 0.5])
                for line in ax.lines
            )
            np.testing.assert_allclose(ax.axis(), [0, T_FINAL, 0, NUM_E + NUM_I + 1])
        assert [ax.get_title() for ax in fig.axes] == [
            rf'$\overline{{I}}_E={drive:g}$' for drive in I_EXT_E_VEC
        ]
        assert fig.axes[-1].get_xlabel() == '$t$ [ms]'
    finally:
        ns.plt.close(fig)
