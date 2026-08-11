import ast
from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "31_ING_Rhythms"


def test_ing_entraining_e_cells_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial E and I
    # firing, most cells active) rather than exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "ING_ENTRAINING_E_CELLS" / "main.py")
    assert len(py.t_e_spikes) > 1000
    assert len(py.t_i_spikes) > 200
    counts_e = np.bincount(py.i_e_spikes.astype(int), minlength=py.num_e)
    counts_i = np.bincount(py.i_i_spikes.astype(int), minlength=py.num_i)
    assert np.count_nonzero(counts_e) > 0.9 * py.num_e
    assert np.count_nonzero(counts_i) > 0.9 * py.num_i


def test_ing_entraining_e_cells_2_structure():
    py = load_python_port(
        ROOT / "python" / PYTHON_BASE / "ING_ENTRAINING_E_CELLS_2" / "main.py"
    )
    assert not hasattr(py, "results")
    results = py.run_drive_panels(t_final_run=100.0)
    assert len(results) == 3
    for t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes in results:
        assert len(t_e_spikes) > 500
        assert len(t_i_spikes) > 200
        assert len(t_e_spikes) == len(i_e_spikes)
        assert len(t_i_spikes) == len(i_i_spikes)


def assert_exact_main_guard(path):
    tree = ast.parse(path.read_text(), filename=str(path))
    guards = [
        node
        for node in tree.body
        if isinstance(node, ast.If)
        and any(
            isinstance(candidate, ast.Name) and candidate.id == "__name__"
            for candidate in ast.walk(node.test)
        )
    ]
    expected = ast.parse('if __name__ == "__main__":\n    main()\n').body[0]
    assert len(guards) == 1
    assert ast.dump(guards[0], include_attributes=False) == ast.dump(
        expected, include_attributes=False
    )


def test_ing_entraining_e_cells_2_default_main_preserves_reference_plot(
    monkeypatch,
):
    path = (
        ROOT
        / "python"
        / PYTHON_BASE
        / "ING_ENTRAINING_E_CELLS_2"
        / "main.py"
    )
    py = load_python_port(path)
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

    monkeypatch.setitem(py.main.__globals__, "run_drive_panels", fake_run_drive_panels)
    saved = []
    monkeypatch.setattr(py.plt, "savefig", saved.append)

    py.main()
    fig = py.plt.gcf()
    try:
        assert calls == [()]
        assert saved == ["fig.png"]
        assert len(fig.axes) == 3
        for ax in fig.axes:
            np.testing.assert_array_equal(
                ax.get_yticks(), [py.num_i, py.num_e + py.num_i]
            )
            assert any(
                np.array_equal(
                    line.get_ydata(), [py.num_i + 0.5, py.num_i + 0.5]
                )
                for line in ax.lines
            )
            np.testing.assert_allclose(
                ax.axis(), [0, py.t_final, 0, py.num_e + py.num_i + 1]
            )
        assert [ax.get_title() for ax in fig.axes] == [
            rf'$\overline{{I}}_E={drive:g}$' for drive in py.i_ext_e_vec
        ]
        assert fig.axes[-1].get_xlabel() == '$t$ [ms]'
    finally:
        py.plt.close(fig)

    assert_exact_main_guard(path)
