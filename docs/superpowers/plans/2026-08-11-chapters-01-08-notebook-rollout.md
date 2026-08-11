# Chapters 01–08 Notebook Rollout Implementation Plan

> **For agentic workers:** Execute chapter-by-chapter, one commit per chapter. This is a lighter-weight plan than the chapter-9 pilot: one review/verification pass per chapter (not per example), no new MATLAB-comparison tests for examples that don't already have one.

**Goal:** Replace the per-example `main.py` scripts in `python/01_*` through `python/08_*` (skipping 02/06, which don't exist) with one `chapterNN.ipynb` per chapter, mirroring the `python/chapter09.ipynb` pilot pattern.

**Architecture:** Each chapter's example folders collapse into one notebook with a `simulate_*`/plot pair per example (or per brian-mirrored group, where `brian/chapterNN.ipynb` already reuses one function across several figures — copy that decomposition instead of a rigid 1-folder-1-function mapping). Existing MATLAB-comparison tests get repointed at the notebook via `load_notebook_definitions_as_module`; examples with no prior test get a cheap smoke test only (call `simulate_*()` with defaults, assert it runs and returns finite arrays) — no new MATLAB reference work in this pass.

**Tech Stack:** Python, Jupyter/nbformat, ipywidgets, pytest, existing `tests/matlab_ref.py` helpers.

## Global Constraints

- Do **not** move non-RTM gating functions into `mnd.core` — `mnd.core.alpha_h/alpha_m/.../n_inf` are already reserved for the RTM family (used by ch09 and by ch05's `RTM_VOLTAGE_TRACE`/`RTM_GATING_VARIABLES`, byte-identical, safe to import). Every other gating family (classical-HH ×3 variants, Wang-Buzsáki ×2, Erisir) stays as a private `def` inside its own notebook — same function names are fine since each notebook is its own namespace.
- Only `python/05_.../RTM_VOLTAGE_TRACE` and `python/05_.../THREE_MODELS_GATING_VARIABLES/RTM_GATING_VARIABLES` import gating functions from `mnd.core`; everything else keeps its own local defs.
- Where a `brian/chapterNN.ipynb` exists, mirror its function decomposition (which examples share one `simulate_*`) and its Colab-install first cell (via `python scripts/add_colab_cells.py python/chapterNN.ipynb`). Chapters 03 and (partially) 08 (`QIF_INFINITE_THRESHOLD`, `THREE_CIRCLES`) have no brian analog — use ch09's one-function-per-example convention for those.
- Keep each chapter's `README.md` in `tests/test_python_chapter_docs.py`'s `CHAPTERS` set (do not exclude, unlike ch09) — rewrite its `## Code examples` / `## Running the examples` sections to reference `chapterNN.ipynb` instead of the deleted subfolders, and keep every required heading intact.
- Delete each chapter's example subfolders (`main.py` + `fig*.png`) only after the notebook and tests are green.
- No new `matlab/` reference work and no new RMSE tolerances — only repoint tests that already exist today.

---

## Task 1: Chapter 01 — `HH_VOLTAGE_TRACE`

**Files:**
- Create: `python/chapter01.ipynb`
- Delete: `python/01_Modeling_a_Single_Neuron/HH_VOLTAGE_TRACE/` (main.py + fig_1_3.png)
- Modify: `python/01_Modeling_a_Single_Neuron/README.md`
- Test: `tests/test_ch01_hh_voltage_trace.py` (new — no prior plain test existed)

**Steps:**
- [ ] Write `chapter01.ipynb`: Colab-install cell, one markdown section, `simulate_hh_voltage_trace(c=1.0, g_k=36.0, g_na=120.0, g_l=0.3, v_k=-82.0, v_na=45.0, v_l=-59.0, i_ext=7.0, t_final=200.0, dt=0.01)` porting `HH_VOLTAGE_TRACE/main.py`'s classical-HH gating + `odeint` run, returning `(t, v)`; `plot_hh_voltage_trace(t, v)`; an `interact(...)` cell.
- [ ] Add `tests/test_ch01_hh_voltage_trace.py`: `ns = load_notebook_definitions_as_module(ROOT/"python"/"chapter01.ipynb")`; `t, v = ns.simulate_hh_voltage_trace()`; assert `len(t) == len(v)` and `np.all(np.isfinite(v))` and `v.max() > 0` (it should spike).
- [ ] Run `uv run pytest tests/test_ch01_hh_voltage_trace.py tests/test_ch01_brian_hh_voltage_trace.py -v` (the brian test is unaffected, just confirming nothing else broke) — fix until green.
- [ ] Update `README.md`: replace `HH_VOLTAGE_TRACE/` references with a link to `chapter01.ipynb`; update `## Running the examples` to `jupyter notebook chapter01.ipynb`.
- [ ] Delete `python/01_Modeling_a_Single_Neuron/HH_VOLTAGE_TRACE/`.
- [ ] Run `uv run pytest tests/test_python_chapter_docs.py -v -k ch01_or_all` (or just the full file) to confirm doc-coverage tests still pass.
- [ ] Commit: `git add python/chapter01.ipynb python/01_Modeling_a_Single_Neuron tests/test_ch01_hh_voltage_trace.py && git commit -m "feat: add python/chapter01.ipynb, replace HH_VOLTAGE_TRACE main.py"`

## Task 2: Chapter 03 — `HH_GATING_VARIABLES`

**Files:**
- Create: `python/chapter03.ipynb`
- Delete: `python/03_The_Classical_HH_ODEs/HH_GATING_VARIABLES/`
- Modify: `python/03_The_Classical_HH_ODEs/README.md`
- Test: `tests/test_ch03_hh_gating_variables.py` (new)

**Steps:**
- [ ] Write `chapter03.ipynb` (no brian mirror to copy — use ch09's single-function convention): `simulate_hh_gating_variables(v_range=np.arange(-100, 50, 0.01))` returning `(v, m_inf_v, h_inf_v, n_inf_v, tau_m, tau_h, tau_n)`, porting the classical-HH gating functions (keep the file's unconditional-division `alpha_m`, no guard — it's fine for array input); `plot_hh_gating_variables(...)` for the 3×2 grid; `interact(...)` cell (sliders less meaningful here since there's no free params beyond the v-range — a simple call is enough, interact optional if there's nothing to tune).
- [ ] Add `tests/test_ch03_hh_gating_variables.py`: load notebook definitions, call `simulate_hh_gating_variables()`, assert all returned `m_inf/h_inf/n_inf` arrays are within `[0, 1]` and `tau_*` arrays are positive and finite.
- [ ] Run `uv run pytest tests/test_ch03_hh_gating_variables.py -v` — fix until green.
- [ ] Update README, delete old folder, rerun `tests/test_python_chapter_docs.py`.
- [ ] Commit.

## Task 3: Chapter 04 — `HH_LIMIT_CYCLE`, `HH_REFRACTORINESS`, `HH_SOLUTION`

**Files:**
- Create: `python/chapter04.ipynb`
- Delete: `python/04_Numerical_Solution_of_HH_ODEs/{HH_LIMIT_CYCLE,HH_REFRACTORINESS,HH_SOLUTION}/`
- Modify: `python/04_Numerical_Solution_of_HH_ODEs/README.md`
- Modify (repoint, keep RMSE tolerances as-is): `tests/test_ch04_hh_limit_cycle.py`, `tests/test_ch04_hh_refractoriness.py`, `tests/test_ch04_hh_solution.py`

**Steps:**
- [ ] Write `chapter04.ipynb` with 3 `simulate_*` functions (keep them separate, matching the 3 existing tests, even though brian's chapter04 reuses one `simulate_HH_neuron` — the 3 python tests each need their own callable):
  - `simulate_hh_limit_cycle(i_ext=10.0, t_final=50.0, dt=0.01)` → `(t, v, n)` (v=-50, m=m_inf(v), h=0.6, n=0.4 initial state, classical-HH gating, guarded `alpha_m`).
  - `simulate_hh_refractoriness(i_ext=10.0, t_final=50.0, dt=0.01, pulse_onset=9.0)` → `(t, v)` (same gating/derivative pattern as `HH_REFRACTORINESS/main.py`, `derivative(x0, t, i_ext, pulse_onset)` signature, single pulse-onset run — default matches what the existing test exercises: `args=(10.0, 9.0)`).
  - `simulate_hh_solution(i_ext=10.0, t_final=50.0, dt=0.01)` → `(t, v, m, h, n)`.
  - Plot functions + `interact()` per example, matching each figure's original layout.
- [ ] Repoint the 3 test files: replace `load_python_port(...)` + manual `odeint` calls with `ns = load_notebook_definitions_as_module(ROOT/"python"/"chapter04.ipynb")` and the matching `ns.simulate_*()` call, keeping the same `trace_rmse` assertions/tolerances (`< 1.0` for v, `< 0.01` for n).
- [ ] Run `uv run pytest tests/test_ch04_hh_limit_cycle.py tests/test_ch04_hh_refractoriness.py tests/test_ch04_hh_solution.py -v -m ""` (MATLAB comparison — run if a MATLAB license is available in this environment; otherwise run without `-m ""` so they skip cleanly, and note in the commit message that MATLAB comparison wasn't re-verified locally).
- [ ] Update README, delete old folders, rerun `test_python_chapter_docs.py`.
- [ ] Commit.

## Task 4: Chapter 05 — RTM / WB / Erisir voltage traces + gating variables

**Files:**
- Create: `python/chapter05.ipynb`
- Delete: `python/05_The_Simple_Model_of_Neurons_in_Rodent_Brains/{RTM_VOLTAGE_TRACE,WB_VOLTAGE_TRACE,ERISIR_VOLTAGE_TRACE,ERISIR_VOLTAGE_TRACE_2,THREE_MODELS_GATING_VARIABLES}/`
- Modify: `python/05_The_Simple_Model_of_Neurons_in_Rodent_Brains/README.md`
- Modify (repoint): `tests/test_ch05_rtm_voltage_trace.py`, `tests/test_ch05_erisir_voltage_trace_2.py`
- Test: `tests/test_ch05_wb_voltage_trace.py`, `tests/test_ch05_wb_voltage_trace_1996.py`, `tests/test_ch05_erisir_voltage_trace.py`, `tests/test_ch05_rtm_gating_variables.py` (new smoke tests — no prior plain tests existed for these)

**Steps:**
- [ ] Write `chapter05.ipynb`:
  - `simulate_rtm_voltage_trace(c=1, g_k=80, g_na=100, g_l=0.1, v_k=-100, v_na=50, v_l=-67, i_ext=1.5, t_final=100, dt=0.01)` → `(t, v)`, importing gating functions **from `mnd.core`** (byte-identical family).
  - `simulate_rtm_gating_variables(...)` mirroring `THREE_MODELS_GATING_VARIABLES/RTM_GATING_VARIABLES/main.py`, also importing from `mnd.core`.
  - `simulate_wb_voltage_trace(c=1.0, g_k=9.0, g_na=35.0, g_l=0.1, v_k=-90.0, v_na=55.0, v_l=-65.0, i_ext=0.75, t_final=100.0, dt=0.01)` → `(t, v)`, private WB gating (2003-style, no φ).
  - `simulate_wb_voltage_trace_1996(..., phi=5.0)` → `(t, v)`, private WB-1996 gating (distinct coefficients from the function above — keep both, do not merge).
  - `simulate_erisir_voltage_trace(c=1, g_k=224.0, g_na=112, g_l=0.5, v_k=-90.0, v_na=60, v_l=-70, i_ext=7.0, t_final=100, dt=0.01, n_power=2)` → `(t, v)`, single function with `n_power` kwarg (mirrors `brian/chapter05.ipynb`'s `simulate_Erisir_neuron`), covering both `ERISIR_VOLTAGE_TRACE` (`n_power=2`) and `ERISIR_VOLTAGE_TRACE_2` (`n_power=4`).
  - Plot + `interact()` per example.
- [ ] Repoint `test_ch05_rtm_voltage_trace.py` → `ns.simulate_rtm_voltage_trace()`, keep `rmse < 2.0`.
- [ ] Repoint `test_ch05_erisir_voltage_trace_2.py` → `ns.simulate_erisir_voltage_trace(n_power=4)`, keep `rmse < 5.0`.
- [ ] Add smoke tests (finite-arrays-and-spikes-look-reasonable, same style as Task 1/2) for `simulate_wb_voltage_trace`, `simulate_wb_voltage_trace_1996`, `simulate_erisir_voltage_trace(n_power=2)`, `simulate_rtm_gating_variables`.
- [ ] Run `uv run pytest tests/test_ch05_*.py -v` (skip `-m ""` unless MATLAB is available) — fix until green.
- [ ] Update README, delete old folders, rerun `test_python_chapter_docs.py`.
- [ ] Commit.

## Task 5: Chapter 07 — LIF + HH-subthreshold examples

**Files:**
- Create: `python/chapter07.ipynb`
- Delete: `python/07_Linear_Integrate_and_Fire_(LIF)_Neurons/{LIF_NEURON_WITH_HH,LIF_VOLTAGE_TRACE,LIF_VOLTAGE_TRACE_2,SUBTHR_FOR_HH,TAU_M_FOR_HH}/`
- Modify: `python/07_Linear_Integrate_and_Fire_(LIF)_Neurons/README.md`
- Test: `tests/test_ch07_hh_subthreshold.py`, `tests/test_ch07_lif_voltage_trace.py`, `tests/test_ch07_lif_voltage_trace_2.py` (new smoke tests — no prior plain tests existed)

**Steps:**
- [ ] Write `chapter07.ipynb`, mirroring `brian/chapter07.ipynb`'s 2-function decomposition:
  - `simulate_hh_subthreshold(i_ext=7.0, t_final=100.0, dt=0.01)` → `(t, v, m, h, n)`, one classical-HH run (private gating, unguarded `alpha_m` as in the source files) that backs all three of `LIF_NEURON_WITH_HH`, `SUBTHR_FOR_HH`, `TAU_M_FOR_HH`'s figures.
  - `plot_lif_neuron_with_hh(t, v)`, `plot_subthr_for_hh(t, v, m, h, n, g_na=120.0, g_k=36.0, g_l=0.3, v_na=45.0, v_k=-82.0, v_l=-59.0)` (derives `I_na/I_k/I_l/I_tot`), `plot_tau_m_for_hh(t, m, h, n, g_na=120.0, g_k=36.0, g_l=0.3)` (derives `tau`) — three plot functions consuming the one `simulate_hh_subthreshold` output, per markdown section.
  - `simulate_lif_neuron(tau_m=10.0, i=0.11, t_final=100.0, dt=0.01)` → `(t, v)`, RK4 stepper + reset loop (port `integrate_rk4`), backing both `LIF_VOLTAGE_TRACE` (`tau_m=10, I=0.11`) and `LIF_VOLTAGE_TRACE_2` (`tau_m=2, I=1/(1-exp(-20/tau_m))/tau_m`) via different kwargs.
  - `interact()` cells for `simulate_hh_subthreshold` and `simulate_lif_neuron`.
- [ ] Add smoke tests: `simulate_hh_subthreshold()` returns finite arrays and spikes; `simulate_lif_neuron(tau_m=10, i=0.11)` and `simulate_lif_neuron(tau_m=2, i=...)` return finite, bounded-in-`[0,1]`-ish (reset at 1) traces.
- [ ] Run `uv run pytest tests/test_ch07_*.py -v` — fix until green.
- [ ] Update README, delete old folders, rerun `test_python_chapter_docs.py`.
- [ ] Commit.

## Task 6: Chapter 08 — QIF / theta examples

**Files:**
- Create: `python/chapter08.ipynb`
- Delete: `python/08_Quadratic_Integrate_and_Fire_(QIF)_and_Theta_Neurons/{QIF_INFINITE_THRESHOLD,QIF_VOLTAGE_TRACE,THETA_FIRING,THREE_CIRCLES}/`
- Modify: `python/08_Quadratic_Integrate_and_Fire_(QIF)_and_Theta_Neurons/README.md`
- Modify (repoint): `tests/test_ch08_qif_infinite_threshold.py`
- Test: `tests/test_ch08_qif_voltage_trace.py`, `tests/test_ch08_theta_firing.py` (new smoke tests)

**Steps:**
- [ ] Write `chapter08.ipynb`:
  - `simulate_qif_infinite_threshold(tau_m=2., i=0.15)` → closed-form `(T, t_ast, v_0_to_1, v_1_to_inf, v_minus_inf_to_0)`, porting the analytic derivation verbatim.
  - `simulate_qif_voltage_trace(tau_m=2., i=0.15, t_final=150, dt=0.01)` → `(t, v)`, Heun/RK2 loop with reset (drop the incremental in-loop plotting and the stray `exit()` — just integrate and return arrays; plotting moves to a separate `plot_qif_voltage_trace(t, v)`).
  - `simulate_theta_firing(tau_m=0.5, i=0.505, t_final=150, dt=0.001)` → `(t, theta)`, `plot_theta_firing(t, theta)` plots `1 - cos(theta)`.
  - `plot_three_circles()` — schematic only, no `simulate_*` needed, ports `THREE_CIRCLES/main.py` using `mnd.core.draw_arrow` (already shared, no duplication).
  - `interact()` cells for `simulate_qif_voltage_trace` and `simulate_theta_firing`.
- [ ] Repoint `test_ch08_qif_infinite_threshold.py` → `ns.simulate_qif_infinite_threshold()`, keep the tight `rtol=1e-6`/`atol=1e-6` closed-form assertions.
- [ ] Add smoke tests for `simulate_qif_voltage_trace` and `simulate_theta_firing` (finite arrays, values land in the expected bounded range).
- [ ] Run `uv run pytest tests/test_ch08_*.py -v -m ""` (or without `-m ""` if no MATLAB) — fix until green.
- [ ] Update README, delete old folders, rerun `test_python_chapter_docs.py`.
- [ ] Commit.

---

## Final check (after all 6 chapters)

- [ ] `uv run pytest tests/ -v -k "ch01 or ch03 or ch04 or ch05 or ch07 or ch08"` full pass.
- [ ] `uv run pytest tests/test_python_chapter_docs.py -v` full pass (doc-coverage gate for all 37 chapters, unaffected by exclusion this time).
- [ ] `uv run ruff check mnd/ tests/`.
- [ ] Manually spot-open 1–2 of the new notebooks in Jupyter to confirm `interact()` widgets render (not required per-chapter, just once at the end, per the "lighter process" request).
