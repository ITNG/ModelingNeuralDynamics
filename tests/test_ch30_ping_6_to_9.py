import ast
from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "30_The_PING_Model_of_Gamma_Rhythms"


def test_ping_6_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing in
    # all three connectivity panels) rather than exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_6" / "main.py")
    for suffix in ("_1", "_2", "_3"):
        assert not hasattr(py, f"t_e_spikes{suffix}")
        assert not hasattr(py, f"t_i_spikes{suffix}")
    panels = py.run_connectivity_panels(t_final_run=200.0)
    assert len(panels) == 3
    for t_e, i_e, t_i, i_i in panels:
        assert len(t_e) > 1000
        assert len(t_i) > 200
        assert len(t_e) == len(i_e)
        assert len(t_i) == len(i_i)


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


def test_ping_6_default_main_preserves_reference_plot(monkeypatch):
    path = ROOT / "python" / PYTHON_BASE / "PING_6" / "main.py"
    py = load_python_port(path)
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
        py.main.__globals__, "run_connectivity_panels", fake_run_connectivity_panels
    )
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
                ax.axis(),
                [py.t_final - 200, py.t_final, 0, py.num_e + py.num_i + 1],
            )
        assert fig.axes[-1].get_xlabel() == '$t$ [ms]'
    finally:
        py.plt.close(fig)

    assert_exact_main_guard(path)


def test_ping_7_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_7" / "main.py")
    assert len(py.t_e_spikes) > 200
    assert len(py.t_i_spikes) > 30
    assert py.lfp.min() > -100 and py.lfp.max() < 50


def test_ping_8_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_8" / "main.py")
    assert len(py.t_e_spikes) > 200
    assert len(py.t_i_spikes) > 30
    assert py.lfp.min() > -100 and py.lfp.max() < 50


def test_ping_9_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_9" / "main.py")
    assert len(py.t_e_spikes) > 200
    assert len(py.t_i_spikes) > 30
    assert py.lfp.min() > -100 and py.lfp.max() < 50
