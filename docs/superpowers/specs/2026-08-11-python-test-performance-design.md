# Python Notebook Test Performance Design

## Goal

Reduce the wall-clock time of tests for notebooks under `python/` while
preserving meaningful numerical coverage. Brian notebooks and Brian-marked
tests are explicitly out of scope.

## Current State

The non-Brian tests do not execute whole Python notebooks. There are 121 test
call sites using `load_notebook_definitions_as_module`, which parses each
notebook and executes only imports plus function and class definitions. The
same notebook is nevertheless parsed, transformed, compiled, and initialized
repeatedly across tests. Some long-running simulation kernels already use
Numba, while many notebook imports include interactive dependencies that tests
do not need.

## Design

### Definitions-only loading

Tests must continue to avoid executing notebook demonstration, plotting,
widget, and other top-level runtime cells. The loader will cache safe,
immutable extraction or compilation work by resolved notebook path and file
metadata. Each request will still receive an isolated namespace unless tests
prove that sharing a namespace is safe.

Imports needed by extracted definitions will remain available. Imports used
only by interactive top-level cells may be excluded when static analysis and
tests show that no extracted function depends on them.

### Benchmark-driven Numba use

The suite will be profiled with the project interpreter at
`/home/ziaee/envs/mnd/bin/python`. Numba will be added only to numerical kernels
whose measured runtime is dominated by compatible Python loops and whose
end-to-end test time improves after including compilation cost. Array-oriented
NumPy/SciPy code and short functions will not be JIT-decorated without evidence
of a net gain.

Notebook-facing APIs and default scientific resolutions will remain unchanged
unless a production-kernel optimization is behaviorally equivalent.

### Smaller test workloads

Individual tests may use shorter simulation durations, smaller parameter
sweeps, or coarser grids when their assertions retain the same behavioral
contract. Full-resolution defaults remain available to notebook users. Tests
that validate precise reference traces or bifurcation boundaries will retain
the resolution required by their tolerances.

## Workflow

1. Establish per-test and suite-level baselines using the project interpreter.
2. Rank non-Brian Python-notebook tests by wall-clock duration.
3. Optimize the definitions loader and verify isolation and extraction rules.
4. Address remaining hotspots one at a time using either a smaller test
   workload or a measured Numba-compatible kernel optimization.
5. Re-run affected numerical assertions after every change.
6. Compare final wall-clock timings with the baseline.

## Verification

- No non-Brian test uses the whole-notebook execution helper for a path under
  `python/`.
- Loader tests prove that top-level runtime cells are not executed and that
  repeated loads avoid repeated parsing/compilation without leaking mutable
  state between namespaces.
- Each optimized hotspot retains its existing numerical assertions.
- The default non-Brian suite passes with no unknown-marker warnings.
- Before/after duration reports demonstrate a net suite improvement; changes
  with no measurable gain are reverted.

## Non-goals

- Optimizing or executing anything under `brian/`.
- Weakening numerical assertions solely to make tests pass.
- Adding Numba decorators broadly without profiling evidence.
- Changing notebook presentation or interactive controls.
