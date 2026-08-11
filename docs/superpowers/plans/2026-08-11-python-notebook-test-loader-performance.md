# Python Notebook Test Loader Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make non-Brian tests load Python notebook definitions faster while guaranteeing that they never execute whole notebooks or share mutable execution namespaces.

**Architecture:** Split notebook-definition loading into a cached compilation stage and an uncached namespace-execution stage. Cache compiled code by resolved path plus file metadata, remove the UI-only `ipywidgets` import from extracted definitions, and keep each returned namespace isolated. After this phase, profile the full non-Brian suite and use the measured duration table to write a second, file-specific plan for Numba and smaller test workloads.

**Tech Stack:** Python 3.11+, pytest, nbformat, `functools.lru_cache`, Python AST, Jupyter notebooks

## Global Constraints

- Only tests for notebooks under `python/` are in scope.
- Anything under `brian/` and every `brian2`-marked test is out of scope.
- Tests must not execute demonstration, plotting, widget, or other top-level notebook runtime cells.
- Each load receives an isolated namespace; only immutable compiled code is cached.
- Use `/home/ziaee/envs/mnd/bin/python` for benchmarks and verification.
- Add Numba only after end-to-end timings show a net gain including compilation cost.
- Preserve full-resolution notebook defaults; smaller workloads belong in tests unless a production-kernel change is behaviorally equivalent.

---

## File Structure

- Modify `tests/matlab_ref.py`: isolate and cache notebook-definition compilation while retaining the public `load_notebook_definitions_as_module(path)` interface.
- Create `tests/test_notebook_definition_loader.py`: cover compilation caching, cache invalidation, namespace isolation, UI-import exclusion, and top-level-cell exclusion.
- No notebook is modified in this phase. Notebook changes belong to the measured hotspot phase.

### Task 1: Cache Definitions-Only Notebook Compilation

**Files:**
- Modify: `tests/matlab_ref.py:1-10,162-198`
- Create: `tests/test_notebook_definition_loader.py`

**Interfaces:**
- Consumes: notebook paths accepted by the existing `load_notebook_definitions_as_module(path)` helper.
- Produces: `_compile_notebook_definitions(path: Path, mtime_ns: int, size: int) -> CodeType`, an internal `lru_cache(maxsize=64)` compilation function.
- Preserves: `load_notebook_definitions_as_module(path) -> SimpleNamespace` with a fresh namespace on every call.

- [ ] **Step 1: Record the uncached repeated-load baseline**

Run:

```bash
PYTHONPATH=tests /home/ziaee/envs/mnd/bin/python - <<'PY'
from pathlib import Path
from time import perf_counter

from matlab_ref import load_notebook_definitions_as_module

path = Path("python/chapter17.ipynb")
started = perf_counter()
namespaces = [load_notebook_definitions_as_module(path) for _ in range(10)]
elapsed = perf_counter() - started
assert len({id(ns) for ns in namespaces}) == 10
print(f"10 repeated definition loads: {elapsed:.6f}s")
PY
```

Expected: ten distinct namespaces are created; retain the printed elapsed time for comparison in Step 7.

- [ ] **Step 2: Write the failing cache and isolation tests**

Create `tests/test_notebook_definition_loader.py`:

```python
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
```

- [ ] **Step 3: Run the cache test to verify RED**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest tests/test_notebook_definition_loader.py -v
```

Expected: `test_repeated_loads_reuse_compilation_but_isolate_namespaces` fails because `nbformat.read` is called twice. The invalidation test should pass as a characterization of current uncached behavior.

- [ ] **Step 4: Implement cached compilation with fresh namespace execution**

In `tests/matlab_ref.py`, add these imports:

```python
from functools import lru_cache
from types import CodeType, SimpleNamespace
```

Replace the compilation portion of `load_notebook_definitions_as_module` with:

```python
@lru_cache(maxsize=64)
def _compile_notebook_definitions(
    path: Path, mtime_ns: int, size: int
) -> CodeType:
    del mtime_ns, size  # Values participate in the cache key.
    import nbformat

    notebook = nbformat.read(path, as_version=4)
    nodes = []
    for cell in notebook.cells:
        if cell["cell_type"] != "code":
            continue
        tree = ast.parse("".join(cell["source"]), filename=str(path))
        nodes.extend(
            node
            for node in tree.body
            if isinstance(
                node,
                (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
            )
            and not (
                isinstance(node, ast.ImportFrom)
                and node.module == "ipywidgets"
            )
            and not (
                isinstance(node, ast.Import)
                and any(alias.name == "ipywidgets" for alias in node.names)
            )
        )
    module = ast.Module(body=nodes, type_ignores=[])
    return compile(ast.fix_missing_locations(module), str(path), "exec")


def load_notebook_definitions_as_module(path):
    """Load imports and definitions without running notebook examples."""
    import matplotlib

    matplotlib.use("Agg")

    path = Path(path).resolve()
    stat = path.stat()
    code = _compile_notebook_definitions(path, stat.st_mtime_ns, stat.st_size)
    namespace = {"__name__": path.stem}
    cwd = os.getcwd()
    try:
        os.chdir(path.parent)
        exec(code, namespace)
    finally:
        os.chdir(cwd)
    return SimpleNamespace(**namespace)
```

Do not cache `SimpleNamespace` or `namespace`; executing the cached `CodeType` creates new function objects bound to a new globals dictionary.

- [ ] **Step 5: Run the cache tests to verify GREEN**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest tests/test_notebook_definition_loader.py -v
```

Expected: both tests pass.

- [ ] **Step 6: Verify representative existing notebook tests**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch01_hh_voltage_trace.py \
  tests/test_ch09_notebook_smoke.py \
  tests/test_ch13_notebook_smoke.py \
  tests/test_ch17_legacy_f_i_curves.py \
  tests/test_ch20_rtm_plot_q.py \
  tests/test_ch22_wilson_cowan_phase_plane.py
```

Expected: all selected non-Brian tests pass. No `ipywidgets` or full-notebook runtime cell is required.

- [ ] **Step 7: Repeat the benchmark and require a net improvement**

Run the exact Step 1 command in a fresh process.

Expected: ten distinct namespaces are still created, and repeated-load time is lower than the Step 1 baseline. If it is not lower across three runs, revert the cache implementation and retain only the regression tests that guard whole-notebook exclusion.

- [ ] **Step 8: Commit the loader optimization**

```bash
git add tests/matlab_ref.py tests/test_notebook_definition_loader.py
git commit -m "perf: cache Python notebook definition compilation"
```

### Task 2: Guard Against Whole Python Notebook Execution

**Files:**
- Modify: `tests/test_notebook_definition_loader.py`

**Interfaces:**
- Consumes: `tests/` Python source files and the loader names defined in `tests/matlab_ref.py`.
- Produces: regression tests that reject top-level notebook execution and Python-notebook use of `load_notebook_as_module`.

- [ ] **Step 1: Add a top-level execution regression test**

Append:

```python
def test_loader_does_not_execute_top_level_runtime_cells(tmp_path):
    path = tmp_path / "chapter.ipynb"
    write_notebook(
        path,
        "def answer():\n    return 42",
        "raise RuntimeError('top-level notebook cell executed')",
    )

    namespace = load_notebook_definitions_as_module(path)

    assert namespace.answer() == 42
```

- [ ] **Step 2: Add a UI-only import exclusion test**

Append:

```python
def test_loader_excludes_ipywidgets_imports(tmp_path):
    path = tmp_path / "chapter.ipynb"
    write_notebook(
        path,
        "from ipywidgets import interact\n\ndef answer():\n    return 42",
    )

    namespace = load_notebook_definitions_as_module(path)

    assert namespace.answer() == 42
    assert not hasattr(namespace, "interact")
```

- [ ] **Step 3: Add a repository-wide static guard for Python notebook tests**

Append:

```python
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
```

This test deliberately ignores Brian tests and catches the established call form used throughout this repository.

- [ ] **Step 4: Run the guard tests**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest tests/test_notebook_definition_loader.py -v
```

Expected: all five loader tests pass.

- [ ] **Step 5: Commit the execution guards**

```bash
git add tests/test_notebook_definition_loader.py
git commit -m "test: prevent full Python notebook execution"
```

### Task 3: Profile the Remaining Non-Brian Suite

**Files:**
- Read: `tests/test_*.py`
- Read: `python/chapter*.ipynb`
- Do not modify production or test files in this task.

**Interfaces:**
- Consumes: the optimized loader from Task 1 and default pytest marker exclusions.
- Produces: a ranked duration table used as the sole input to the follow-up Numba/workload implementation plan.

- [ ] **Step 1: Verify collection excludes Brian tests**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest --collect-only -q tests/test_*_brian_*.py
```

Expected: zero selected tests and all Brian cases reported as deselected.

- [ ] **Step 2: Run the non-Brian suite with durations**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q --durations=50
```

Expected: the suite passes and prints the 50 slowest setup/call/teardown durations. Do not add Numba based on source appearance alone.

- [ ] **Step 3: Re-run the ten slowest call phases individually**

For each node ID printed as a `call` duration in Step 2, run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q '<exact-node-id>' --durations=1
```

Expected: retain only hotspots that reproduce within 15% across three runs. Classify each retained hotspot as notebook loading, simulation loop, SciPy solver, or oversized test workload.

- [ ] **Step 4: Create the measured follow-up plan**

Use `superpowers:writing-plans` to create a second plan whose tasks name the exact slow test, notebook function, baseline duration, proposed optimization, numerical assertion, and required after-time. Use smaller test workloads before Numba when both preserve the same contract. Do not include Brian files.

The follow-up plan is required because the exact notebook files eligible for Numba cannot be chosen responsibly until Steps 2-3 provide timings.
