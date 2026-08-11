import ast
import os
from pathlib import Path

import nbformat
import pytest

from matlab_ref import load_notebook_definitions_as_module


def whole_notebook_loader_references(source: str) -> list[int]:
    """Return lines that reference the whole-notebook execution helper."""
    loader_name = "load_notebook_" + "as_module"
    references = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.ImportFrom) and any(
            alias.name == loader_name for alias in node.names
        ):
            references.add(node.lineno)
        elif isinstance(node, ast.Name) and node.id == loader_name:
            references.add(node.lineno)
        elif isinstance(node, ast.Attribute) and node.attr == loader_name:
            references.add(node.lineno)
    return sorted(references)


def non_brian_whole_notebook_loader_violations(tests_dir: Path) -> list[str]:
    violations = []
    for test_path in sorted(tests_dir.glob("test_*.py")):
        if "_brian_" in test_path.name:
            continue
        if whole_notebook_loader_references(test_path.read_text()):
            violations.append(test_path.name)
    return violations


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
    assert first.answer.__globals__ is not second.answer.__globals__
    first.answer.__globals__["namespace_sentinel"] = object()
    assert "namespace_sentinel" not in second.answer.__globals__


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


@pytest.mark.parametrize(
    "source",
    [
        """\
from matlab_ref import load_notebook_as_module
namespace = load_notebook_as_module (
    ROOT / 'python' / 'chapter.ipynb'
)
""",
        """\
from matlab_ref import load_notebook_as_module as execute_notebook
namespace = execute_notebook(ROOT / 'python' / 'chapter.ipynb')
""",
        """\
loader = load_notebook_as_module
namespace = loader(ROOT / 'python' / 'chapter.ipynb')
""",
        """\
import matlab_ref as helpers
namespace = helpers.load_notebook_as_module(ROOT / 'python' / 'chapter.ipynb')
""",
    ],
    ids=["spacing-multiline-single-quotes", "import-alias", "variable", "module-alias"],
)
def test_whole_notebook_loader_guard_catches_bypass_forms(source):
    assert whole_notebook_loader_references(source)


def test_whole_notebook_loader_guard_preserves_brian_exclusion(tmp_path):
    source = """\
from matlab_ref import load_notebook_as_module
namespace = load_notebook_as_module(ROOT / 'python' / 'chapter.ipynb')
"""
    (tmp_path / "test_ch01_python.py").write_text(source)
    (tmp_path / "test_ch01_brian_python.py").write_text(source)

    assert non_brian_whole_notebook_loader_violations(tmp_path) == [
        "test_ch01_python.py"
    ]


def test_python_notebook_tests_never_use_whole_notebook_loader():
    tests_dir = Path(__file__).resolve().parent
    assert non_brian_whole_notebook_loader_violations(tests_dir) == []
