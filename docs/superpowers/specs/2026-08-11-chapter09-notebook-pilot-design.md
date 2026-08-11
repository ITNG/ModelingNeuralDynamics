# Chapter 9 notebook pilot — design

## Context

`ModelingNeuralDynamics` currently has two ported tracks of the book's MATLAB
originals: `python/` (plain NumPy/SciPy scripts, one self-contained
`main.py` per example, grouped into per-chapter folders) and `brian/` (one
notebook per chapter, some chapters already Colab-ready).

The goal is to make the `python/` track runnable and explorable on Colab too,
with `ipywidgets` sliders so a reader can play with model parameters
in-browser. Converting all 39 `python/` chapters (254 `main.py` files, ~150+
dependent test files) in one pass is too large to design or execute safely at
once. This spec covers a single pilot chapter — **Chapter 9, Spike Frequency
Adaptation** — to work out the notebook structure, the shared-code strategy,
and the test-adaptation pattern that later chapters will reuse.

Chapter 9 was chosen because: it already has a Brian2 notebook
(`brian/chapter09.ipynb`) to mirror conventions from, it has real
duplication (four examples share nine HH gating functions verbatim) to
validate the shared-module decision against, and its guide doc
(`python/09_Spike_Frequency_Adaptation/README.md`) is already excluded from
`tests/test_python_chapter_docs.py`'s coverage checks (deferred by an earlier
effort), so this pilot doesn't collide with that in-flight work.

## Scope

**In scope:** chapter 9 only — replace its 8 `main.py` scripts with one
notebook, update the 6 existing test files that exercise them, add the
shared gating functions to `mnd`, add `ipywidgets` as an explicit dependency,
touch up chapter 9's README.

**Out of scope:** any other chapter's `python/` scripts, the `brian/` track,
the deferred chapter-9 guide-doc effort (`python/README.md` index,
`test_python_chapter_docs.py` coverage), and the CI-runtime work from the
current session (already merged separately). The chapter-by-chapter rollout
plan for the remaining 38 chapters is a follow-up, scoped once this pilot is
validated.

## Current state (chapter 9)

Eight example folders under `python/09_Spike_Frequency_Adaptation/`, each
with a standalone `main.py`:

| Folder | What it computes |
|---|---|
| `ADAPTATION_MAP` | spike-to-spike adaptation map φ(z) |
| `CALCIUM_RISE` | voltage-dependent calcium target c∞(v) |
| `LIF_ADAPT` | LIF neuron with reset-based adaptation |
| `M_CURRENT` | M-gate steady state and time constant |
| `RTM_AHP` | RTM neuron with calcium-dependent AHP current |
| `RTM_AHP_RESTING` | same, at zero external drive |
| `RTM_M` | RTM neuron with M-current |
| `RTM_M_RESTING` | same, at zero external drive |
| `V_V_TILDE` | subthreshold v vs. its adaptation-free counterpart |

`RTM_AHP`, `RTM_AHP_RESTING`, `RTM_M`, `RTM_M_RESTING` each redefine the same
nine HH gating functions (`alpha_h`, `alpha_m`, `alpha_n`, `beta_h`,
`beta_m`, `beta_n`, `h_inf`, `m_inf`, `n_inf`) byte-for-byte.

Six test files (`tests/test_ch09_{rtm_m,rtm_ahp,adaptation_map,
calcium_rise,lif_adapt,m_current,v_v_tilde}.py`) import each `main.py` via
`load_python_port()` and check the Python output against MATLAB references
(skipped in CI without a MATLAB license, still exercised on release tags per
the earlier CI change).

## Design

### 1. Notebook layout

New file: `python/chapter09.ipynb`, flat under `python/` (mirrors
`brian/chapter09.ipynb`'s placement directly under `brian/`, not nested in a
chapter subfolder). The 8 example folders and their `main.py`/`fig.png`
files are deleted — the notebook is the sole source for chapter 9 going
forward.

First cell: the same Colab-install pattern already used in
`brian/chapter09.ipynb`:

```python
import subprocess, sys
if "google.colab" in sys.modules:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "modelingneuraldynamics"])
```

One markdown section per example (mirroring the book's structure), each
containing:
- a `simulate_*` function taking the example's tunable parameters as keyword
  arguments with the book's original defaults, returning the arrays needed
  for plotting and testing
- a `plot_*` function (or inline plotting) for the static figure
- an `interact(...)` cell wiring `ipywidgets` sliders to the `simulate_*` +
  `plot_*` pair for interactive exploration

This is the same one-function-per-example convention
`brian/chapter09.ipynb` already uses (e.g. `simulate_RTM_M_neuron`), chosen
over flat top-level variables (the `main.py` scripts' current style) for two
reasons: consistency with the existing Brian2 notebook, and because
`ipywidgets.interact()` needs a function to attach sliders to regardless, so
keeping parameters as function arguments avoids a separate wrapper layer.

Example signature:

```python
def simulate_rtm_m(i_ext=1.5, g_m=0.25, t_final=300, dt=0.01):
    """RTM neuron with an M-current."""
    ...
    return t, v, w
```

### 2. Shared code → `mnd`

The nine gating functions move into the `mnd` package (alongside the
existing plain-NumPy helpers in `mnd/core.py`), since they're duplicated
byte-for-byte across four of chapter 9's own examples today. `mnd` is
already a real installable package (pip-installable as
`modelingneuraldynamics`, used by `brian/chapter09.ipynb`'s Colab-install
cell), so this also solves Colab portability: a notebook opened directly in
Colab only fetches the single `.ipynb` file, not sibling files, so shared
code must come from something `pip install`-able rather than a colocated
`lib.py`.

`chapter09.ipynb` imports these from `mnd` instead of redefining them per
example. Exact submodule placement (e.g. `mnd.core` vs. a new
`mnd.core.hh_gating`) is an implementation detail decided while writing the
code, matching whatever the existing `mnd/core.py` organization suggests.

### 3. `ipywidgets` dependency

Add `ipywidgets` explicitly to `pyproject.toml`'s dependencies. It's
currently only present transitively (via `jupyter`); the notebook code
imports it directly (`from ipywidgets import interact`), so it should be a
direct, explicit dependency.

### 4. Test adaptation

All 6 test files switch from:

```python
py = load_python_port(ROOT / "python" / PYTHON_DIR / "main.py")
```

to:

```python
ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
```

`load_notebook_definitions_as_module` already exists in `tests/matlab_ref.py`
for exactly this case: it parses each code cell with `ast`, keeps only
`import`/`def`/`class` nodes, and executes just those — so the notebook's
`interact(...)` cells (top-level statements) are never executed during
tests. This is why `load_notebook_definitions_as_module` is used here rather
than `load_notebook_as_module` (which the existing `brian/` chapter 9 tests
use, full-notebook exec) — chapter 9's python notebook deliberately adds
widget cells that `load_notebook_as_module` would otherwise run.

Each test then calls the relevant `ns.simulate_*(...)` function directly
with default arguments (matching the book's parameters) instead of
re-running `odeint` against a module-level `derivative`/`x0` pair. The
MATLAB-comparison logic itself (RMSE tolerances, spike-time comparisons)
does not change — only how the Python-side trace is obtained.

Example (`test_ch09_rtm_m.py`), before/after shape:

```python
# before
py = load_python_port(ROOT / "python" / "09_Spike_Frequency_Adaptation" / "RTM_M" / "main.py")
t_p = py.np.arange(0, py.t_final, py.dt)
sol = odeint(py.derivative, py.x0, t_p)
v_p, w_p = sol[:, 0], sol[:, 3]

# after
ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
t_p, v_p, w_p = ns.simulate_rtm_m()
```

### 5. Docs

`python/09_Spike_Frequency_Adaptation/README.md` is already excluded from
`test_python_chapter_docs.py`'s `CHAPTERS` coverage set (deferred by the
earlier guide-docs effort), so no test depends on its content matching the
current folder layout. It will still be updated (not left stale) to
reference the new notebook instead of the deleted per-example folders and
`python main.py` instructions, since it remains a real, linkable file in the
repo even though unvalidated.

## Testing

- Run `uv run pytest tests/test_ch09_*.py -v -m ""` locally (the pilot
  machine has a MATLAB license) to confirm the RMSE/spike-time checks still
  pass against MATLAB with the new notebook-based trace extraction.
- Run `uv run ruff check .` (notebook code isn't linted directly, but the
  `mnd` gating-function move is plain `.py`).
- Manually open `python/chapter09.ipynb` in Jupyter to confirm the
  `interact()` widgets render and each example's static figure still
  matches the deleted `main.py`'s `fig.png` output visually.
- Confirm `load_notebook_definitions_as_module` does *not* trigger any
  `interact()` call during test collection (i.e. tests run without a
  Jupyter frontend/comm target present).

## Follow-up (not in this pilot)

- Rollout plan for the remaining 38 chapters, once this pattern is
  validated.
- Whether/how the deferred chapter-9 guide doc
  (`python/09_Spike_Frequency_Adaptation/README.md` → linked from
  `python/README.md`) gets un-deferred now that the notebook exists.
