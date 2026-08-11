# Final Review Fix Report

Date: 2026-08-11

Branch: `perf/python-notebook-tests`

Worktree: `/home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics/.worktrees/perf-python-notebook-tests`

## Outcome

All requested final-review corrections were made without touching Brian files or implementing any planned production optimization:

- Replaced the exact-text whole-notebook regression guard with an AST identifier/import-level prohibition for every non-Brian `test_*.py` file.
- Added direct bypass coverage for whitespace plus multiline/single-quote calls, import aliases, variable indirection, and module aliases.
- Preserved and directly tested the `_brian_` filename exclusion.
- Strengthened namespace isolation with an exact `__globals__` identity assertion and a mutation-leak assertion.
- Rewrote follow-up Task 4 so both legacy scripts produce `main()`, their `__main__` guards call it exactly, and their existing plotting blocks move intact.
- Restored the PING y-ticks and the ING population separator plus y-ticks in the implementation snippets.
- Added exact implementation-ready tests that mock only the expensive runner, call the real default `main()` path with no shortened-duration argument, inspect real Matplotlib axes, and structurally verify the `__main__` guard.

## Root-Cause Verification

The original guard used this single lexical shape:

```python
if "load_notebook_as_module(ROOT / \"python\"" in source:
```

It therefore depended on an exact function spelling immediately followed by `(`, exact spacing, and double quotes. AST inspection confirmed that aliases and variable references were semantic uses that the substring could not see.

`tests/matlab_ref.py::load_python_port` calls `spec_from_file_location(path.stem, path)`. Both target files are named `main.py`, so their imported `__name__` is `main`, not `__main__`. The former Task 4 verification could pass without entering either guard. Reading the current plotting blocks also confirmed the omitted behavior:

- `PING_6/main.py`: `ax.set_yticks([num_i, num_e + num_i])`
- `ING_ENTRAINING_E_CELLS_2/main.py`: the population separator line and the same y-tick setup

## Changed Files

- `tests/test_notebook_definition_loader.py`
  - AST reference detector and non-Brian file scanner
  - four bypass-form parameter cases
  - Brian-exclusion fixture
  - globals identity and mutation-leak assertions
- `docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md`
  - corrected Task 4 interfaces, plotting snippets, mocked default-main tests, structural main-guard assertions, commands, and expected results
- `.superpowers/sdd/2026-08-11-python-notebook-test-loader-performance/final-fix-report.md`
  - this evidence report

No path under `brian/` and no Brian test file changed.

## TDD Evidence

### Behavior-preserving extraction before RED

The original substring was first extracted into `whole_notebook_loader_references` and the repository scanner into `non_brian_whole_notebook_loader_violations` without changing behavior.

Command:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_notebook_definition_loader.py::test_python_notebook_tests_never_use_whole_notebook_loader
```

Exact output:

```text
.                                                                        [100%]
1 passed in 1.67s
```

### RED: bypass forms and Brian-exclusion fixture

After adding the desired cases but before replacing the substring detector:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_notebook_definition_loader.py \
  -k 'bypass_forms or preserves_brian_exclusion'
```

Exact pytest result:

```text
FFFFF                                                                    [100%]
=================================== FAILURES ===================================
_ test_whole_notebook_loader_guard_catches_bypass_forms[spacing-multiline-single-quotes] _
E       assert []
_____ test_whole_notebook_loader_guard_catches_bypass_forms[import-alias] ______
E       assert []
_______ test_whole_notebook_loader_guard_catches_bypass_forms[variable] ________
E       assert []
_____ test_whole_notebook_loader_guard_catches_bypass_forms[module-alias] ______
E       assert []
__________ test_whole_notebook_loader_guard_preserves_brian_exclusion __________
E       AssertionError: assert [] == ['test_ch01_python.py']
=========================== short test summary info ============================
FAILED tests/test_notebook_definition_loader.py::test_whole_notebook_loader_guard_catches_bypass_forms[spacing-multiline-single-quotes]
FAILED tests/test_notebook_definition_loader.py::test_whole_notebook_loader_guard_catches_bypass_forms[import-alias]
FAILED tests/test_notebook_definition_loader.py::test_whole_notebook_loader_guard_catches_bypass_forms[variable]
FAILED tests/test_notebook_definition_loader.py::test_whole_notebook_loader_guard_catches_bypass_forms[module-alias]
FAILED tests/test_notebook_definition_loader.py::test_whole_notebook_loader_guard_preserves_brian_exclusion
5 failed, 5 deselected in 2.68s
```

The failures were the expected missing semantic-detection behavior, not collection or syntax errors.

### GREEN: AST guard

After implementing the minimal AST detector for `ast.ImportFrom`, `ast.Name`, and `ast.Attribute` references:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_notebook_definition_loader.py
```

Exact output:

```text
..........                                                               [100%]
10 passed in 2.70s
```

### Namespace-isolation mutation check

Namespace isolation already existed, so no production fix was required. To prove the new assertion catches the specific leak, `load_notebook_definitions_as_module` was temporarily mutated to execute repeated loads in one shared globals dictionary while still returning distinct `SimpleNamespace` objects and distinct function objects.

Command:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_notebook_definition_loader.py::test_repeated_loads_reuse_compilation_but_isolate_namespaces
```

Exact salient output:

```text
F                                                                        [100%]
=================================== FAILURES ===================================
_________ test_repeated_loads_reuse_compilation_but_isolate_namespaces _________
>       assert first.answer.__globals__ is not second.answer.__globals__
E       assert {...} is not {...}
=========================== short test summary info ============================
FAILED tests/test_notebook_definition_loader.py::test_repeated_loads_reuse_compilation_but_isolate_namespaces
1 failed in 1.79s
```

The failure occurred on the new `__globals__` identity assertion after the pre-existing `first is not second` and `first.answer is not second.answer` assertions passed. The deliberate mutation was then fully restored; `git diff --exit-code -- tests/matlab_ref.py` emitted no output and exited 0.

### Final fresh GREEN

Command:

```bash
git diff --exit-code -- tests/matlab_ref.py && \
/home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_notebook_definition_loader.py && \
git status --short && \
git diff --check
```

Exact output:

```text
..........                                                               [100%]
10 passed in 2.71s
```

The two Git diff checks and `git status --short` emitted no output at that point.

## Plan Validations

### Python-fence and Task 4 interface checks

Command:

```bash
/home/ziaee/envs/mnd/bin/python - <<'PY'
from pathlib import Path
import re

path = Path('docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md')
text = path.read_text()
task4 = text.split('### Task 4:', 1)[1].split('### Task 5:', 1)[0]
blocks = re.findall(r'```python\n(.*?)\n```', task4, re.S)
for index, block in enumerate(blocks, 1):
    compile(block, f'{path}:task4-block-{index}', 'exec')
print(f'Task 4 Python fences compile: {len(blocks)}/{len(blocks)}')
assert task4.count('ax.set_yticks([num_i, num_e + num_i])') == 2
assert task4.count("ax.plot([0, t_final], [num_i + 0.5, num_i + 0.5], '--k', linewidth=1)") == 2
assert task4.count('def assert_exact_main_guard(path):') == 2
assert task4.count('assert calls == [()]') == 2
assert task4.count('\n    py.main()\n') == 2
assert task4.count('if __name__ == "__main__":\n    main()') == 2
assert '`load_python_port` imports with `__name__ == "main"`' in task4
assert 'Move the existing `fig, axes` block into `main()` intact' in task4
assert 'Move the existing plotting block into `main()` intact' in task4
print('Task 4 default-main/plot interface checks: PASS')
PY
```

Exact output:

```text
Task 4 Python fences compile: 6/6
Task 4 default-main/plot interface checks: PASS
```

### Markdown fence balance

Command:

```bash
/home/ziaee/envs/mnd/bin/python - <<'PY'
from pathlib import Path
path = Path('docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md')
count = sum(line.startswith('```') for line in path.read_text().splitlines())
print(f'{path}: fence lines={count}; balanced={count % 2 == 0}')
assert count % 2 == 0
PY
```

Exact output:

```text
docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md: fence lines=62; balanced=True
```

### Placeholder scan

Command:

```bash
rg -n 'TBD|TODO|implement later|fill in details|Similar to|\.\.\.|# existing|first_kernel|threshold corrections|simulate_from_rest|rest_duration|driven_duration' \
  docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md
```

Exact output: none; exit code 1, meaning no placeholder pattern matched. The combined validation wrapper printed:

```text
placeholder scan: PASS (no matches)
```

### Diff and scope checks

Commands:

```bash
git diff --check
git diff --name-only HEAD
```

Exact output before the plan commit:

```text
docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md
```

`git diff --check` emitted no output and exited 0. At the final pre-report test point, `git status --short` and `git diff --check` both emitted no output.

## Commits

- `418c893 test: harden whole-notebook loader guard`
- `7721400 docs: verify default plotting entry points`

This report is committed separately after its contents are finalized; the resulting report commit is returned to the parent agent together with the two fix commits.

## Self-Review

- **Finding 1:** The scanner now reasons over parsed identifiers/imports, not quote or whitespace spelling. Every requested bypass form failed before the AST change and passed afterward. The repository-level guard parses all current non-Brian tests successfully.
- **Finding 2:** Task 4 explicitly documents why `load_python_port` cannot verify a main guard. Each script now has a planned `main()` interface, a no-argument runner call, an exact AST guard assertion, and a lightweight real-plot test with only the expensive simulation and file save replaced.
- **Plot preservation:** PING contains y-ticks, separator, 200ms window, and three-panel plotting. ING contains separator, y-ticks, full default axis, drive titles, and three-panel plotting. Both save `fig.png`.
- **Namespace isolation:** Distinct function globals and absence of a sentinel injected into the first globals dictionary are asserted. A shared-globals mutation fails exactly at the new identity assertion.
- **Scope:** No Brian file and no production optimization file changed. `tests/matlab_ref.py` was changed only for the deliberate mutation check and restored byte-for-byte before final verification.
- **Plan quality:** All six Task 4 Python fences compile; the placeholder scan is empty; markdown fences are balanced; planned producer/consumer names and tuple orders agree.

## Concerns

No blocking concerns. Task 4 remains an unexecuted follow-up plan by design; this fix wave corrected and validated the plan but did not implement its PING/ING production refactors or run the long bounded simulations. The default-main tests specified by the plan intentionally use real Matplotlib objects and fake only the expensive runner plus `savefig`, so their future execution will validate plotting behavior without producing published files.
