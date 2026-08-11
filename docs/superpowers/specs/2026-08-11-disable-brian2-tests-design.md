# Temporarily Exclude Brian2 Tests by Default

## Goal

Make the normal pytest run faster while Brian2 development is paused, without
deleting or permanently disabling any Brian2 coverage.

## Design

Pytest collection will assign a `brian2` marker to tests collected from files
whose names match the repository's existing `test_*_brian_*.py` convention. The
marker will be registered in `pyproject.toml`, and the default marker expression
will change from `not slow` to `not slow and not brian2`.

This keeps the current default behavior for slow MATLAB tests and additionally
excludes Brian2 tests. Developers can still run the Brian2 suite explicitly with
`pytest -m brian2`. They can run slow Brian2 tests, if any are added later, with
an explicit expression such as `pytest -m "brian2 and slow"`.

## Scope

- Add one collection hook in `tests/conftest.py` that marks Brian2 test files by
  filename.
- Register and exclude the `brian2` marker by default in `pyproject.toml`.
- Do not edit the 57 individual Brian2 test files.
- Do not alter dependencies or Brian2 production code.

## Verification

- Collection with the default options selects no tests from
  `test_*_brian_*.py` files.
- `pytest -m brian2 --collect-only` selects the Brian2 tests.
- A representative non-Brian2 test still runs normally.
- The pytest configuration emits no unknown-marker warnings.

## Reversal

Remove `and not brian2` from the default marker expression to restore Brian2
tests to the normal suite. The marker registration and collection hook may stay
because they remain useful for targeted runs.
