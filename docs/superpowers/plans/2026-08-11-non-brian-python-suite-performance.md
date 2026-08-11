# Non-Brian Python Suite Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce the measured non-Brian Python test hotspots without weakening their numerical contracts or executing complete notebooks.

**Architecture:** First reduce test-only time horizons and population sizes where measured trials prove the same contract still holds. For legacy scripts, move expensive top-level simulations behind callable entry points so tests can request bounded workloads while normal script execution retains reference defaults. Use cached Numba kernels only for the chapter-35/36 dense f-I scans, where reducing the current grid would remove the plateau-resolution assertion.

**Tech Stack:** Python 3, pytest, NumPy, nbformat notebook definitions, Numba (`njit(cache=True)`), Matplotlib's non-interactive test backend.

## Global Constraints

- Exclude every Brian file and every `brian`-marked test.
- Keep `tests/matlab_ref.py::load_notebook_definitions_as_module` as the only loader used by notebook tests; never call a whole-notebook executor.
- Do not change notebook or script defaults used for published figures.
- Prefer a smaller test workload when it preserves the asserted model behavior; use Numba only for dense scans whose point count is part of the contract.
- Preserve deterministic seeds (`63806`) and existing exact/qualitative numerical assertions.
- Measure call durations with `/home/ziaee/envs/mnd/bin/python -m pytest -q <node> --durations=<N>` on the same machine used for the 2026-08-11 baseline.
- Required after-times are cold-process pytest call durations, not averages from an already-warm interpreter.

## Measured Hotspot Disposition

The retained set is every completed call at or above 5.00s, every cap-active exact call, and the three chapter-17 calls required by the brief's stability reproduction. Stable reproductions use their median once. Completed timings are sorted numerically; cap-active lower bounds precede them because each exceeds the largest completed call. Every retained call has either an implementation task or an explicit deferral below.

| Order | Call | Baseline | Disposition |
|---:|---|---:|---|
| 1 | ch30 `test_ping_6_structure` | >115s exact cap | Task 4; measured 200ms workload preserves all assertions. |
| 2 | ch31 `test_ing_entraining_e_cells_2_structure` | >115s exact cap | Task 4; measured 100ms workload preserves all assertions. |
| 3 | ch35 `test_rtm_f_i_curve_with_inhibition_structure` | >115s exact cap | Task 6; retain the full current grid and accelerate the loop. |
| 4 | ch35 `test_rtm_f_i_curve_with_inhibition_2_structure` | >115s exact cap | Task 6; retain the full current grid and accelerate the loop. |
| 5 | ch36 `test_rtm_f_i_curve_pulsed_excitation_2_structure` | >115s exact cap | Task 6; retain the full plateau-resolution grid and accelerate the loop. |
| 6 | ch36 `test_rtm_f_i_curve_pulsed_excitation_structure` | >113s contextual lower bound | Task 6; retain the full plateau-resolution grid and accelerate the loop. |
| 7 | ch24 `test_rtm_e_to_e_heterogeneous_is_sane` | 107.86s | Task 3; measured `N=4,t_final=100` workload preserves the contract. |
| 8 | ch30 `test_ping_5_structure` | 94.17s | Defer: one run and no bounded-workload trial; profile after cap-active PING 6 with its seeded contract isolated. |
| 9 | ch32 `test_m_current_ping_1_from_rest_structure` | 78.30s | Defer: no shortened from-rest trial; its two-phase history needs independent measurements before changing the horizon. |
| 10 | ch32 `test_m_current_ping_1_structure` | 74.18s | Task 5; measured 200ms workload preserves activity and voltage assertions. |
| 11 | ch32 `test_m_current_ping_3_closeup_structure` | 72.46s | Defer: no shortened trial; the PING-1 measurement is not evidence for this variant. |
| 12 | ch35 `test_periodic_inhibition_f_i_curve_structure` | 64.08s | Task 6; retain the full grid and accelerate the loop. |
| 13 | ch31 `test_ing_entraining_e_cells_structure` | 42.33s | Defer: one run and no shortened-workload trial; E2 is prioritized because it exceeded the cap and has measured bounds. |
| 14 | ch32 `test_poisson_ping_2_structure` | 38.28s | Defer: below the 45s network-task ceiling and no shortened-workload evidence. |
| 15 | ch38 `test_poisson_ping_population_statistics_match_matlab` | 36.71s | Defer: its statistical reference contract requires a separate reproducibility study before reducing samples. |
| 16 | ch32 `test_poisson_ping_3_structure` | 35.34s | Defer: below 45s and no shortened-workload evidence. |
| 17 | ch32 `test_poisson_ping_3_voltage_trace_structure` | 34.75s | Defer: below 45s and no shortened-workload evidence. |
| 18 | ch32 `test_ping_clusters_structure` | 32.20s | Defer: below 45s and no shortened-workload evidence. |
| 19 | ch37 `test_ping_thr_1_zoom_structure` | 31.77s | Defer: below 45s and no measured shorter trace that preserves the zoom-window contract. |
| 20 | ch39/40 `test_three_cell_ping_5_structure` | 29.41s | Defer: below 30s and no shortened-workload evidence. |
| 21 | ch33 `test_pinb_1_structure` | 28.24s | Defer: below 30s and no shortened-workload evidence. |
| 22 | ch14 cycle-distance test | 27.06s median | Task 1; measured 200ms workload preserves every assertion. |
| 23 | ch33 `test_pinb_3_structure` | 26.68s | Defer: below 30s and no shortened-workload evidence. |
| 24 | ch33 `test_m_current_ping_7_structure` | 25.92s | Defer: below 30s and no shortened-workload evidence. |
| 25 | ch33 `test_m_current_ping_6_structure` | 25.84s | Defer: below 30s and no shortened-workload evidence. |
| 26 | ch33 `test_m_current_beta_with_gj_structure` | 25.81s | Defer: below 30s and no shortened-workload evidence. |
| 27 | ch33 `test_pinb_2_structure` | 24.55s | Defer: below 30s and no shortened-workload evidence. |
| 28 | ch30 `test_ping_7_structure` | 17.73s | Defer: below 30s and no shortened-workload evidence. |
| 29 | ch30 `test_ping_8_structure` | 16.37s | Defer: below 30s and no shortened-workload evidence. |
| 30 | ch30 `test_ping_9_structure` | 14.99s | Defer: below 30s and no shortened-workload evidence. |
| 31 | ch31 `test_ing_10_structure` | 14.88s | Defer: below 30s and no shortened-workload evidence. |
| 32 | ch31 `test_ing_9_structure` | 14.20s | Defer: below 30s and no shortened-workload evidence. |
| 33 | ch31 `test_ing_6_structure` | 13.88s | Defer: below 30s and no shortened-workload evidence. |
| 34 | ch31 `test_ing_5_structure` | 13.87s | Defer: below 30s and no shortened-workload evidence. |
| 35 | ch31 `test_ing_2_structure` | 13.74s | Defer: below 30s and no shortened-workload evidence. |
| 36 | ch31 `test_ing_8_structure` | 13.64s | Defer: below 30s and no shortened-workload evidence. |
| 37 | ch31 `test_ing_4_structure` | 13.58s | Defer: below 30s and no shortened-workload evidence. |
| 38 | ch31 `test_ing_7_structure` | 13.39s | Defer: below 30s and no shortened-workload evidence. |
| 39 | ch31 `test_ing_3_structure` | 13.14s | Defer: below 30s and no shortened-workload evidence. |
| 40 | ch37 `test_ping_thr_1_structure` | 12.99s | Defer: below 30s and no shortened-workload evidence. |
| 41 | ch33 `test_m_current_ping_5_structure` | 11.62s | Defer: below 30s and no shortened-workload evidence. |
| 42 | ch33 `test_m_current_ping_4_structure` | 11.29s | Defer: below 30s and no shortened-workload evidence. |
| 43 | ch24 `test_rtm_e_to_e_network_smoke` | 9.98s | Defer: already below 10s. |
| 44 | ch24 `test_rtm_splay_smoke` | 8.71s | Defer: already below 10s. |
| 45 | ch38 `test_poisson_ping_is_reproducible_and_active` | 7.41s | Defer: already below 10s. |
| 46 | ch38 phase-response mismatch case | 6.78s | Defer: already below 10s. |
| 47 | ch18 `test_h_current_matches_matlab` | 6.58s | Defer: already below 10s and exact MATLAB-reference comparison. |
| 48 | ch38 phase-response plus case | 6.58s | Defer: already below 10s. |
| 49 | ch18 `test_plot_modified_tau_r_matches_matlab` | 6.55s | Defer: already below 10s and exact MATLAB-reference comparison. |
| 50 | ch23 `test_plot_f_entrainment_matches_matlab` | 6.38s | Defer: already below 10s and exact MATLAB-reference comparison. |
| 51 | ch09 notebook smoke | 6.24s | Defer: already below 10s. |
| 52 | ch17 HH legacy f-I | 5.35s median | Task 2; measured 100ms workload preserves the qualitative response. |
| 53 | ch17 RTM onset legacy f-I | 1.11s median | Task 2; retained because the brief required its three-run reproduction. |
| 54 | ch17 RTM legacy f-I | 0.79s median | Task 2; retained because the brief required its three-run reproduction. |

Additional context-only cap-active calls were not completed and therefore cannot be numerically ranked: chapter-32 M-current PING 1 close-up, M-current PING 2 close-up, and Poisson PING 1; chapter-33 M-current PING 8; chapter-35 periodic-inhibition f-I curve 2; chapter-37 thresholding; the final chapter-38 phase-statistics case; and chapter-40 PING with STDP. They are explicitly deferred until each can be isolated under the 120-second profiling cap; no implementation target is asserted from incomplete evidence.

---

### Task 1: Bound the Chapter-14 Cycle-Distance Smoke Workload

**Measured baseline and target:** `tests/test_ch14_hh_reduced_cycle_distance.py::test_hh_reduced_cycle_distance_fixed_points_are_sane` reproduced at 27.59s, 26.55s, and 27.06s (27.06s median). A measured `t_final=200.0` trial retained all three fixed-point and window assertions and took 6.64s wall time. Required pytest call time: **< 8.0s**.

**Files:**
- Modify: `tests/test_ch14_hh_reduced_cycle_distance.py:22-33`
- Read only: `python/chapter14.ipynb` functions `simulate_hh_reduced_cycle_distance`, `hh_reduced_simulate`, and `hh_reduced_find_fixed_point`

**Interfaces:**
- Consumes: `simulate_hh_reduced_cycle_distance(i_ext_vec=(5.5, 5.4, 5.3), t_final=1000.0, dt=0.01) -> list[tuple]`
- Produces: the same three panels and assertions, using an explicit test-only `t_final=200.0`

- [ ] **Step 1: Change the test to state its bounded workload explicitly**

```python
def test_hh_reduced_cycle_distance_fixed_points_are_sane():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter14.ipynb")
    panels = ns.simulate_hh_reduced_cycle_distance(t_final=200.0)

    assert [panel[0] for panel in panels] == [5.5, 5.4, 5.3]
    for i_ext, v_c, n_c, v_attr, n_attr, v_rep, n_rep in panels:
        assert -68 < v_c < -65
        assert 0.35 < n_c < 0.38
        assert _in_window(v_attr, n_attr).any()
        assert _in_window(v_rep, n_rep).any()
```

- [ ] **Step 2: Run the exact test and enforce the after-time**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch14_hh_reduced_cycle_distance.py::test_hh_reduced_cycle_distance_fixed_points_are_sane \
  --durations=1
```

Expected: `1 passed`; call duration below 8.0s; three currents remain `[5.5, 5.4, 5.3]`; every attracting and repelling trajectory enters the plotted window.

- [ ] **Step 3: Commit the isolated workload change**

```bash
git add tests/test_ch14_hh_reduced_cycle_distance.py
git commit -m "test: bound chapter 14 cycle-distance workload"
```

---

### Task 2: Bound and Strengthen the Chapter-17 Legacy f-I Smoke Tests

**Measured baselines and targets:** Across three runs, HH was 5.35s/5.20s/5.52s (5.35s median), RTM onset was 1.11s/1.09s/1.15s (1.11s median), and RTM was 0.79s/0.78s/0.83s (0.79s median). A measured combined `t_final=100.0` trial took 3.99s and returned HH `[0, 62.9167, 75.1984]`, HH backward `[0, 62.9162, 74.6513]`, RTM `[0, 0, 42.6533]`, RTM backward `[0, 0, 42.6524]`, and three zero onset frequencies. Required file time: **< 5.0s**; HH call **< 3.5s**; each RTM call **< 0.8s**.

**Files:**
- Modify: `tests/test_ch17_legacy_f_i_curves.py:15-39`
- Read only: `python/chapter17.ipynb` functions `_legacy_scan`, `simulate_hh_f_i_curve_legacy`, `simulate_rtm_f_i_curve_legacy`, and `simulate_rtm_f_i_curve_at_onset_legacy`

**Interfaces:**
- Consumes: the three existing notebook simulation functions and their current tuple return types
- Produces: bounded 100ms scans plus qualitative firing assertions, rather than shape/finite checks alone

- [ ] **Step 1: Bound HH and assert its silent-to-firing transition**

```python
f_forward, f_backward, i_ext_vec = ns.simulate_hh_f_i_curve_legacy(
    i_ext_vec=np.array([3.0, 8.0, 13.0]), t_final=100.0
)
assert f_forward.shape == f_backward.shape == i_ext_vec.shape == (3,)
assert np.all(np.isfinite(f_forward)) and np.all(np.isfinite(f_backward))
assert f_forward[0] == f_backward[0] == 0.0
assert f_forward[-1] > 70.0 and f_backward[-1] > 70.0
```

- [ ] **Step 2: Bound RTM and assert its high-current firing response**

```python
f_forward, f_backward, i_ext_vec = ns.simulate_rtm_f_i_curve_legacy(
    i_ext_vec=np.array([0.0, 0.5, 1.0]), t_final=100.0
)
assert f_forward.shape == f_backward.shape == i_ext_vec.shape == (3,)
assert np.all(np.isfinite(f_forward)) and np.all(np.isfinite(f_backward))
assert f_forward[0] == f_backward[0] == 0.0
assert f_forward[-1] > 40.0 and f_backward[-1] > 40.0
```

- [ ] **Step 3: Bound onset and assert the expected short-window result**

```python
f_forward, i_ext_vec, I_c, C = ns.simulate_rtm_f_i_curve_at_onset_legacy(
    n_points=3, t_final=100.0
)
assert f_forward.shape == i_ext_vec.shape == (3,)
assert np.array_equal(f_forward, np.zeros(3))
assert I_c is None and C is None
```

- [ ] **Step 4: Run all three tests and enforce the after-times**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch17_legacy_f_i_curves.py --durations=3
```

Expected: `3 passed in < 5.0s`; HH call below 3.5s; both RTM calls below 0.8s.

- [ ] **Step 5: Commit the chapter-17 workload change**

```bash
git add tests/test_ch17_legacy_f_i_curves.py
git commit -m "test: bound chapter 17 legacy f-i scans"
```

---

### Task 3: Use a Representative Chapter-24 Heterogeneous Network

**Measured baseline and target:** `tests/test_ch24_rtm_e_to_e_heterogeneous.py::test_rtm_e_to_e_heterogeneous_is_sane` took 107.86s with `N=30, t_final=200.0`. A measured `N=4, t_final=100.0` trial took 19.07s, produced two spikes from every neuron, preserved the seeded drive/coupling ranges, and preserved zero diagonal coupling. Required pytest call time: **< 22.0s**.

**Files:**
- Modify: `tests/test_ch24_rtm_e_to_e_heterogeneous.py:10-28`
- Read only: `python/chapter24.ipynb` function `simulate_rtm_e_to_e_heterogeneous`

**Interfaces:**
- Consumes: `simulate_rtm_e_to_e_heterogeneous(N=30, t_final=200.0, dt=0.01, seed=63806)`; all omitted model parameters retain their notebook defaults
- Produces: a deterministic four-neuron, 100ms smoke workload with the same shape/range/connectivity contract and two spikes per neuron

- [ ] **Step 1: Write the bounded workload and exact numerical assertions**

```python
t_spikes, i_spikes, i_ext, g_syn = ns.simulate_rtm_e_to_e_heterogeneous(
    N=4, t_final=100.0, seed=63806
)
N = len(i_ext)

assert N == 4
assert i_ext.shape == (N,)
assert np.all((0.25 <= i_ext) & (i_ext <= 0.35))
assert g_syn.shape == (N, N)
assert np.all(np.diag(g_syn) == 0)
off_diag = g_syn[~np.eye(N, dtype=bool)]
assert np.all((0.00625 <= off_diag) & (off_diag <= 0.00875))
counts = np.array([np.sum(i_spikes == i) for i in range(1, N + 1)])
assert np.array_equal(counts, np.array([2, 2, 2, 2]))
```

- [ ] **Step 2: Run the exact test and enforce the after-time**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch24_rtm_e_to_e_heterogeneous.py::test_rtm_e_to_e_heterogeneous_is_sane \
  --durations=1
```

Expected: `1 passed`; call duration below 22.0s; seeded spike counts exactly `[2, 2, 2, 2]`.

- [ ] **Step 3: Commit the chapter-24 workload change**

```bash
git add tests/test_ch24_rtm_e_to_e_heterogeneous.py
git commit -m "test: use representative heterogeneous network"
```

---

### Task 4: Separate Reference Runs from PING-6 and ING-Entraining Smoke Runs

**Measured baselines and targets:** PING 6 exceeded the 115s command cap. A definition-only trial of all three existing connectivity panels at 200ms took 38.63s and produced `(E, I)` spike counts `(1888, 437)`, `(1947, 395)`, and `(1943, 424)`, satisfying the existing `>1000`/`>200` assertions. ING entraining E cells 2 exceeded 115s; a definition-only three-drive trial at 100ms took 29.57s and produced `(2059, 510)`, `(2231, 419)`, and `(2341, 447)`, satisfying the existing `>500`/`>200` assertions. Required call times: **PING 6 < 45s** and **ING entraining E cells 2 < 35s**.

**Files:**
- Modify: `python/30_The_PING_Model_of_Gamma_Rhythms/PING_6/main.py`
- Modify: `tests/test_ch30_ping_6_to_9.py:9-18`
- Modify: `python/31_ING_Rhythms/ING_ENTRAINING_E_CELLS_2/main.py`
- Modify: `tests/test_ch31_ing_entraining_e_cells.py:21-26`

**Interfaces:**
- Produces in PING 6: `run_connectivity_panels(t_final_run: float = t_final) -> list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]`
- Produces in ING: `run_drive_panels(t_final_run: float = t_final) -> list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]`
- Normal `python main.py` execution still calls each entry point with the original 2000ms/500ms defaults before plotting.

- [ ] **Step 1: Add the PING-6 callable and move reference execution under the main guard**

```python
def run_connectivity_panels(t_final_run=t_final):
    panels = []
    g1 = make_random_connectivity(p_ee, p_ei, p_ie, p_ii)
    panels.append(simulate(*g1, t_final=t_final_run))

    sparse = (1 / num_e, 1 / num_e, 1 / num_i, 1 / num_i)
    g2 = make_random_connectivity(*sparse)
    panels.append(simulate(*g2, t_final=t_final_run))

    g3 = make_fixed_degree_connectivity(sparse[1], sparse[2], sparse[3], ni=1)
    panels.append(simulate(*g3, t_final=t_final_run))
    return panels


if __name__ == "__main__":
    panels_raw = run_connectivity_panels()
    panels = [(ti, ii, te, ie) for te, ie, ti, ii in panels_raw]
    fig, axes = plt.subplots(3, 1, figsize=(8, 9))
    for ax, (ti, ii, te, ie) in zip(axes, panels):
        if len(ti) > 0:
            ax.plot(ti, ii, '.b', markersize=2)
        if len(te) > 0:
            ax.plot(te, ie + num_i, '.r', markersize=2)
        ax.plot([0, t_final], [num_i + 0.5, num_i + 0.5], '--k', linewidth=1)
        ax.axis([t_final - 200, t_final, 0, num_e + num_i + 1])
    axes[-1].set_xlabel('$t$ [ms]')
    plt.tight_layout()
    plt.savefig("fig.png")
```

Remove the three unconditional panel simulation assignments so importing definitions is cheap and the figure path retains the 2000ms default.

The build/simulate pairs must remain sequential exactly as shown. The module-level RNG is consumed by both connectivity construction and simulation, so constructing every matrix before any simulation would change the seeded reference streams and invalidate the measured spike counts.

- [ ] **Step 2: Update the PING-6 test to call the 200ms workload**

```python
py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_6" / "main.py")
panels = py.run_connectivity_panels(t_final_run=200.0)
assert len(panels) == 3
for t_e, i_e, t_i, i_i in panels:
    assert len(t_e) > 1000
    assert len(t_i) > 200
    assert len(t_e) == len(i_e)
    assert len(t_i) == len(i_i)
```

- [ ] **Step 3: Add the ING callable and move reference execution under the main guard**

```python
def run_drive_panels(t_final_run=t_final):
    return [
        simulate(
            drive * np.ones(num_e) * (1 + sigma_e * rng.standard_normal(num_e)),
            t_final=t_final_run,
        )
        for drive in i_ext_e_vec
    ]


if __name__ == "__main__":
    results = run_drive_panels()
    fig, axes = plt.subplots(3, 1, figsize=(8, 9))
    for ax, drive, (te, ie, ti, ii) in zip(axes, i_ext_e_vec, results):
        if len(ti) > 0:
            ax.plot(ti, ii, '.b', markersize=2)
        if len(te) > 0:
            ax.plot(te, ie + num_i, '.r', markersize=2)
        ax.axis([0, t_final, 0, num_e + num_i + 1])
        ax.set_title(rf'$\overline{{I}}_E={drive:g}$')
    axes[-1].set_xlabel('$t$ [ms]')
    plt.tight_layout()
    plt.savefig("fig.png")
```

Remove the unconditional `results` comprehension so import does not run the 500ms simulations.

- [ ] **Step 4: Update the ING test to call the 100ms workload**

```python
py = load_python_port(
    ROOT / "python" / PYTHON_BASE / "ING_ENTRAINING_E_CELLS_2" / "main.py"
)
results = py.run_drive_panels(t_final_run=100.0)
assert len(results) == 3
for t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes in results:
    assert len(t_e_spikes) > 500
    assert len(t_i_spikes) > 200
    assert len(t_e_spikes) == len(i_e_spikes)
    assert len(t_i_spikes) == len(i_i_spikes)
```

- [ ] **Step 5: Verify both bounded tests and default script smoke behavior**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch30_ping_6_to_9.py::test_ping_6_structure \
  tests/test_ch31_ing_entraining_e_cells.py::test_ing_entraining_e_cells_2_structure \
  --durations=2
```

Expected: `2 passed`; PING 6 below 45s; ING below 35s; the original assertions remain satisfied.

- [ ] **Step 6: Commit the legacy entry-point refactor**

```bash
git add \
  'python/30_The_PING_Model_of_Gamma_Rhythms/PING_6/main.py' \
  tests/test_ch30_ping_6_to_9.py \
  'python/31_ING_Rhythms/ING_ENTRAINING_E_CELLS_2/main.py' \
  tests/test_ch31_ing_entraining_e_cells.py
git commit -m "test: bound legacy network smoke simulations"
```

---

### Task 5: Bound the Measured Chapter-32 M-Current PING-1 Workload

**Measured baseline and target:** `test_m_current_ping_1_structure` took 74.18s. Its measured 200ms trial took 21.81s, produced 368 E spikes and 282 I spikes, and reached `48.377mV`. Required exact-call time: **< 25s**. The from-rest and other close-up variants are excluded because no bounded trial was measured for them; their explicit deferrals are in the disposition table.

**Files:**
- Modify: `python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_1/main.py`
- Modify: `tests/test_ch32_m_current_ping_poisson_ping.py:25-31`

**Interface:**
- Produce exactly `run_smoke(t_final_run: float = t_final) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]` in tuple order `(t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, v_plot)`.
- Normal `python main.py` execution calls `run_smoke()` with the unchanged published default before computing rates and plotting.

- [ ] **Step 1: Add the callable and move the reference run under the main guard**

```python
def run_smoke(t_final_run=t_final):
    connectivity = make_random_connectivity(p_ee, p_ei, p_ie, p_ii)
    return simulate(*connectivity, t_final=t_final_run)


if __name__ == "__main__":
    t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, v_plot = run_smoke()
    f_hat_e = round(len(t_e_spikes) / num_e / t_final * 1000)
    f_hat_i = round(len(t_i_spikes) / num_i / t_final * 1000)
    print(f"f_hat_e = {f_hat_e}")
    print(f"f_hat_i = {f_hat_i}")
```

Move all plotting that consumes these arrays into the same main guard. Loading the script for a test must define functions and constants without constructing connectivity or running `simulate`.

- [ ] **Step 2: Call the measured 200ms workload and preserve the exact observed contract**

```python
py = load_python_port(
    ROOT / "python" / PYTHON_BASE / "M_CURRENT_PING_1" / "main.py"
)
t_e, i_e, t_i, i_i, v_plot = py.run_smoke(t_final_run=200.0)
assert len(t_e) == 368
assert len(t_i) == 282
assert len(t_e) == len(i_e)
assert len(t_i) == len(i_i)
assert np.isclose(v_plot.max(), 48.377, rtol=0.0, atol=0.01)
```

- [ ] **Step 3: Run the exact test and enforce the measured target**

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_1_structure \
  --durations=1
```

Expected: `1 passed`; call below 25s; exact seeded spike totals and voltage maximum remain within the stated tolerance.

- [ ] **Step 4: Commit only the measured variant**

```bash
git add \
  python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_1/main.py \
  tests/test_ch32_m_current_ping_poisson_ping.py
git commit -m "test: bound measured m-current ping workload"
```

---

### Task 6: Make Dense Chapter-35/36 Scans Explicit and Validate Cold JIT Gain

**Measured baselines and targets:** Chapter-35 periodic-inhibition f-I took 64.08s; both RTM inhibition variants exceeded 115s. Chapter-36 pulsed excitation consumed more than 113s for variant 1 and exceeded 115s for variant 2. These five tests require their full 100/101/201-point grids. The companion chapter-35 periodic-inhibition f-I variant 2 remains deferred because profiling produced no attributable duration or lower bound.

| Task 6 target | Measured evidence | Required fresh-cache call time |
|---|---:|---:|
| `test_periodic_inhibition_f_i_curve_structure` | 64.08s completed | <30s and <64.08s |
| `test_rtm_f_i_curve_with_inhibition_structure` | >115s exact cap | <30s and <115s |
| `test_rtm_f_i_curve_with_inhibition_2_structure` | >115s exact cap | <30s and <115s |
| `test_rtm_f_i_curve_pulsed_excitation_structure` | >113s contextual lower bound | <30s and <113s |
| `test_rtm_f_i_curve_pulsed_excitation_2_structure` | >115s exact cap | <30s and <115s |

**Files:**
- Modify: `python/35_Periodic_Inhibition/PERIODIC_INHIBITION_F_I_CURVE/main.py`
- Modify: `python/35_Periodic_Inhibition/RTM_F_I_CURVE_WITH_INHIBITION/main.py`
- Modify: `python/35_Periodic_Inhibition/RTM_F_I_CURVE_WITH_INHIBITION_2/main.py`
- Modify: `python/36_F_I_Curves_Pulsed_Excitation/RTM_F_I_CURVE_PULSED_EXCITATION/main.py`
- Modify: `python/36_F_I_Curves_Pulsed_Excitation/RTM_F_I_CURVE_PULSED_EXCITATION_2/main.py`
- Modify: `tests/test_ch35_periodic_inhibition.py:43-79`
- Modify: `tests/test_ch36_f_i_curves_pulsed_excitation.py:17-38`

**Interfaces:**
- The LIF script produces `compute_f_i_curve(i_values: np.ndarray = I_vec, use_numba: bool = True) -> np.ndarray`.
- RTM inhibition scripts produce `compute_f_i_curves(i_ext_values: np.ndarray = i_ext_vec, use_numba: bool = True) -> tuple[np.ndarray, np.ndarray]` in `(f_vec_tonic, f_vec_periodic)` order.
- Pulsed-excitation scripts produce `compute_f_i_curves(i_ext_values: np.ndarray = i_ext_vec, use_numba: bool = True) -> tuple[np.ndarray, np.ndarray]` in `(f_vec_constant, f_vec_pulsed)` order.
- `load_python_port` executes a module body, so no scan, connectivity construction, or curve array may be created at import time. Full-grid calls and plotting exist only inside `if __name__ == "__main__"`.

- [ ] **Step 1: First make every full scan an explicit callable**

For `PERIODIC_INHIBITION_F_I_CURVE`, move the current top-level nested loop unchanged into `_compute_f_i_curve_python(i_values)`. For the four RTM scripts, rename the existing numerical functions with `_python` suffixes, parameterize their current grids, and have `compute_f_i_curves` call them. Preserve current iteration order and state carry-over.

In `PERIODIC_INHIBITION_F_I_CURVE`, put `f_vec = compute_f_i_curve()` and every existing plot statement beneath its `if __name__ == "__main__"` guard. In each RTM script, put `f_vec_tonic, f_vec_periodic = compute_f_i_curves()` or `f_vec_constant, f_vec_pulsed = compute_f_i_curves()`, as appropriate, plus every existing plot statement beneath that guard.

Update each test to load definitions, verify that legacy result globals are absent, then call the new interface explicitly:

```python
py = load_python_port(path)
assert not hasattr(py, "f_vec")
f_vec = py.compute_f_i_curve()

py = load_python_port(path)
assert not hasattr(py, "f_vec_tonic")
assert not hasattr(py, "f_vec_periodic")
f_vec_tonic, f_vec_periodic = py.compute_f_i_curves()
```

Use corresponding `f_vec_constant`/`f_vec_pulsed` absence checks for chapter 36. Apply the existing full-grid numerical assertions to returned arrays, never import-created globals.

- [ ] **Step 2: Create independent Python and compiled execution paths**

Add these exact imports to each of the five target scripts:

```python
from numba import njit
from numba.extending import register_jitable
```

Do not rebind any Python function name to a Numba dispatcher. Add `@register_jitable` directly above the existing scalar helper definitions (`g`; or `shape`, `m_inf`, `alpha_h`, `beta_h`, `alpha_n`, `beta_n`, `g_periodic`, and `step` as present in that script). This leaves each helper callable by CPython and makes it legal inside compiled outer loops. Do not use a bare `numba` module name anywhere; `njit` and `register_jitable` are the two imported identifiers.

Create dispatchers under new names only: `_compute_f_i_curve_jit = njit(cache=True)(_compute_f_i_curve_python)`, `_f_i_curve_tonic_jit = njit(cache=True)(_f_i_curve_tonic_python)`, `_f_i_curve_periodic_jit = njit(cache=True)(_f_i_curve_periodic_python)`, `_f_i_curve_constant_jit = njit(cache=True)(_f_i_curve_constant_python)`, and `_f_i_curve_pulsed_jit = njit(cache=True)(_f_i_curve_pulsed_python)`, using only the pairs present in each file.

In the LIF public wrapper, choose `_compute_f_i_curve_jit` when `use_numba` is true and `_compute_f_i_curve_python` otherwise, then call the selected kernel with `i_values`. In each RTM public wrapper, choose the two `_jit` names when `use_numba` is true and the two `_python` names otherwise, then call both selected kernels with `i_ext_values` (and the existing periodic amplitude where required). The `use_numba=False` path therefore resolves only the original Python functions and Python-visible registered helpers.

- [ ] **Step 3: Prove Python/JIT equivalence through the public callable**

For each of the five target scripts, compare a three-current probe through the two independent paths:

```python
probe = np.array([grid[0], grid[len(grid) // 2], grid[-1]])
expected = py.compute_f_i_curves(i_ext_values=probe, use_numba=False)
actual = py.compute_f_i_curves(i_ext_values=probe, use_numba=True)
for expected_curve, actual_curve in zip(expected, actual):
    np.testing.assert_allclose(actual_curve, expected_curve, rtol=0.0, atol=1e-10)
```

For `PERIODIC_INHIBITION_F_I_CURVE`, compare `compute_f_i_curve(i_values=probe, use_numba=False)` directly with `compute_f_i_curve(i_values=probe, use_numba=True)`. Keep all existing full-grid assertions: LIF plateau set including `0,40,80,120,160`; RTM 101-point silent-to-active endpoints; pulsed 201-point endpoints and more than ten rounded 40Hz samples.

- [ ] **Step 4: Measure a truly cold cache, then warm reuse**

Create one fresh cache directory and reuse that exact directory for the second command:

```bash
cache_dir=$(mktemp -d /tmp/mnd-numba-task6.XXXXXX)
NUMBA_CACHE_DIR="$cache_dir" /home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch35_periodic_inhibition.py::test_periodic_inhibition_f_i_curve_structure \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_structure \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_2_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_2_structure \
  --durations=5

NUMBA_CACHE_DIR="$cache_dir" /home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch35_periodic_inhibition.py::test_periodic_inhibition_f_i_curve_structure \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_structure \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_2_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_2_structure \
  --durations=5
```

Expected on both runs: `5 passed`; each call satisfies its exact row in the target/evidence/gate table above. Record the generated `cache_dir` path with the timings.

- [ ] **Step 5: Enforce the cold-cost decision gate**

Evaluate each script independently. If its fresh-cache call is not both below 30s and strictly below its measured baseline/lower bound, remove that script's Numba import, registered decorators, dispatcher, and `use_numba` branch. Retain only its import-clean explicit callable and original pure-Python equations. Do not report a Numba speedup for that script. Re-run its numerical test after removal. This gate includes first-call compilation and forbids accepting a warm-only gain.

- [ ] **Step 6: Commit the validated callable/JIT set**

```bash
git add \
  python/35_Periodic_Inhibition/PERIODIC_INHIBITION_F_I_CURVE/main.py \
  python/35_Periodic_Inhibition/RTM_F_I_CURVE_WITH_INHIBITION/main.py \
  python/35_Periodic_Inhibition/RTM_F_I_CURVE_WITH_INHIBITION_2/main.py \
  python/36_F_I_Curves_Pulsed_Excitation/RTM_F_I_CURVE_PULSED_EXCITATION/main.py \
  python/36_F_I_Curves_Pulsed_Excitation/RTM_F_I_CURVE_PULSED_EXCITATION_2/main.py \
  tests/test_ch35_periodic_inhibition.py \
  tests/test_ch36_f_i_curves_pulsed_excitation.py
git commit -m "perf: validate cold dense f-i kernels"
```

---

### Task 7: Verify the Ranked Hotspot Set and the Non-Brian Suite

**Files:**
- Test only: all files changed in Tasks 1-6
- Do not include: `tests/test_*_brian_*.py` or any `brian`-marked node

**Interfaces:**
- Consumes: the bounded workloads and cached kernels from Tasks 1-6
- Produces: a fresh duration table proving every planned target and a bounded-suite result

- [ ] **Step 1: Verify Brian deselection remains intact**

```bash
/home/ziaee/envs/mnd/bin/python -m pytest --collect-only -q tests/test_*_brian_*.py
```

Expected: zero selected, 95 deselected, pytest exit code 5 because no tests were selected.

- [ ] **Step 2: Run the exact optimized hotspot set**

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch14_hh_reduced_cycle_distance.py \
  tests/test_ch17_legacy_f_i_curves.py \
  tests/test_ch24_rtm_e_to_e_heterogeneous.py \
  tests/test_ch30_ping_6_to_9.py::test_ping_6_structure \
  tests/test_ch31_ing_entraining_e_cells.py::test_ing_entraining_e_cells_2_structure \
  tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_1_structure \
  tests/test_ch35_periodic_inhibition.py::test_periodic_inhibition_f_i_curve_structure \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_structure \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_2_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_2_structure \
  --durations=30
```

Expected: all selected tests pass; no call exceeds 45s; every task-specific after-time above is met.

- [ ] **Step 3: Run the remaining non-Brian suite in chapter batches capped at 115 seconds**

Use the same chapter/file batches recorded in `task-3-report.md`. Every command must use `timeout --signal=INT 115s`; if a batch reaches the cap, run its active exact node separately and record it before proceeding to later nodes. The known unrelated documentation inventory failure must be reported separately rather than hidden.

- [ ] **Step 4: Check the final implementation scope and commit history**

```bash
git status --short
git diff --check
git log -7 --oneline
```

Expected: only files named in Tasks 1-6 are changed; no Brian file appears; `git diff --check` emits no output. Do not alter any numerical threshold beyond the exact assertions already specified by Tasks 1-6.
