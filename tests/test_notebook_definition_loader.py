import os
from pathlib import Path

import nbformat

from matlab_ref import load_notebook_definitions_as_module


def write_notebook(path: Path, *sources: str) -> None:
    notebook = nbformat.v4.new_notebook(
        cells=[nbformat.v4.new_code_cell(source) for source in sources]
    )
    nbformat.write(notebook, path)


def test_repeated_loads_reuse_compilation_but_isolate_namespaces(tmp_path, monkeypatch):
    path = tmp_path / "chapter.ipynb"
    write_notebook(path, "def answer():\n    return 42")

    real_read = nbformat.read
    read_count = 0

    def counting_read(*args, **kwargs):
        nonlocal read_count
        read_count += 1
        return real_read(*args, **kwargs)

    monkeypatch.setattr(nbformat, "read", counting_read)
    first = load_notebook_definitions_as_module(path)
    second = load_notebook_definitions_as_module(path)

    assert read_count == 1
    assert first.answer() == second.answer() == 42
    assert first is not second
    assert first.answer is not second.answer


def test_cache_invalidates_when_notebook_changes(tmp_path):
    path = tmp_path / "chapter.ipynb"
    write_notebook(path, "def answer():\n    return 1")
    first = load_notebook_definitions_as_module(path)

    write_notebook(path, "def answer():\n    return 222")
    stat = path.stat()
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1))
    second = load_notebook_definitions_as_module(path)

    assert first.answer() == 1
    assert second.answer() == 222


def test_loader_does_not_execute_top_level_runtime_cells(tmp_path):
    path = tmp_path / "chapter.ipynb"
    write_notebook(
        path,
        "def answer():\n    return 42",
        "raise RuntimeError('top-level notebook cell executed')",
    )

    namespace = load_notebook_definitions_as_module(path)

    assert namespace.answer() == 42


def test_loader_excludes_ipywidgets_imports(tmp_path):
    path = tmp_path / "chapter.ipynb"
    write_notebook(
        path,
        "from ipywidgets import interact\n\ndef answer():\n    return 42",
    )

    namespace = load_notebook_definitions_as_module(path)

    assert namespace.answer() == 42
    assert not hasattr(namespace, "interact")


def test_python_notebook_tests_never_use_whole_notebook_loader():
    tests_dir = Path(__file__).resolve().parent
    violations = []
    for test_path in sorted(tests_dir.glob("test_*.py")):
        if "_brian_" in test_path.name:
            continue
        source = test_path.read_text()
        if "load_notebook_as_module(ROOT / \"python\"" in source:
            violations.append(test_path.name)

    assert violations == []
