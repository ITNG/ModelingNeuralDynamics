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
    g1 = make_random_connectivity(p_ee, p_ei, p_ie, p_ii)
    sparse = (1 / num_e, 1 / num_e, 1 / num_i, 1 / num_i)
    g2 = make_random_connectivity(*sparse)
    g3 = make_fixed_degree_connectivity(sparse[1], sparse[2], sparse[3], ni=1)
    return [simulate(*g, t_final=t_final_run) for g in (g1, g2, g3)]


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

### Task 5: Bound the Chapter-32 M-Current PING Smoke Horizons

**Measured baselines and targets:** M-current PING 1 took 74.18s, from-rest took 78.30s, and PING 3 close-up took 72.46s; PING 1 close-up and PING 2 close-up were cap-active. A measured M-current PING 1 trial at 200ms took 21.81s, produced 368 E spikes and 282 I spikes, and reached `48.38mV`. Required call time for each of the five M-current tests: **< 30s**; combined file segment **< 115s**.

**Files:**
- Modify: `python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_1/main.py`
- Modify: `python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_1_CLOSEUP/main.py`
- Modify: `python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_1_FROM_REST/main.py`
- Modify: `python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_2_CLOSEUP/main.py`
- Modify: `python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_3_CLOSEUP/main.py`
- Modify: `tests/test_ch32_m_current_ping_poisson_ping.py:25-56`

**Interfaces:**
- For each script with a `simulate` function, produce `run_smoke(t_final_run: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]` and keep the current default run inside `if __name__ == "__main__"`.
- For `M_CURRENT_PING_1_FROM_REST`, extract its top-level integration into `simulate_from_rest(t_final_run=t_final, drive_onset=200.0)` returning spike arrays and the tracked voltage.

- [ ] **Step 1: Add callable entry points without changing published defaults**

Use this exact wrapper in the four scripts that already define `simulate`:

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

For the from-rest script, extract the two phases into the following interface. The body calls the already-defined `step` with the shown zero-drive and driven inputs, and uses the same threshold-crossing interpolation as the other four scripts:

```python
def simulate_from_rest(rest_duration=200.0, driven_duration=t_final):
    connectivity = make_random_connectivity(p_ee, p_ei, p_ie, p_ii)
    state = (
        -70.0 * np.ones(num_e),
        m_e_inf(-70.0 * np.ones(num_e)),
        h_e_inf(-70.0 * np.ones(num_e)),
        n_e_inf(-70.0 * np.ones(num_e)),
        np.zeros(num_e), np.zeros(num_e), np.zeros(num_e),
        -75.0 * np.ones(num_i),
        m_i_inf(-75.0 * np.ones(num_i)),
        h_i_inf(-75.0 * np.ones(num_i)),
        n_i_inf(-75.0 * np.ones(num_i)),
        np.zeros(num_i), np.zeros(num_i),
    )
    for _ in range(round(rest_duration / dt)):
        state = step(*state, i_ext_e_rest, i_ext_i_rest, *connectivity)

    driven_e = 3.0 * np.ones(num_e) * (1 + sigma_e * rng.standard_normal(num_e))
    driven_i = 0.7 * np.ones(num_i) * (1 + sigma_i * rng.standard_normal(num_i))
    t_e, i_e, t_i, i_i = [], [], [], []
    for k in range(1, round(driven_duration / dt) + 1):
        old_e, old_i = state[0], state[7]
        state = step(*state, driven_e, driven_i, *connectivity)
        new_e, new_i = state[0], state[7]
        which_e = np.where((old_e > -20) & (new_e <= -20))[0]
        which_i = np.where((old_i > -20) & (new_i <= -20))[0]
        i_e.extend(which_e.tolist())
        i_i.extend(which_i.tolist())
        t_e.extend((((-20 - new_e[which_e]) * (k - 1) * dt
                     + (old_e[which_e] + 20) * k * dt)
                    / (-new_e[which_e] + old_e[which_e])).tolist())
        t_i.extend((((-20 - new_i[which_i]) * (k - 1) * dt
                     + (old_i[which_i] + 20) * k * dt)
                    / (-new_i[which_i] + old_i[which_i])).tolist())
    return tuple(np.asarray(x) for x in (t_e, i_e, t_i, i_i))
```

Normal script execution calls `simulate_from_rest()` with the original 200ms rest phase and 1000ms driven phase.

- [ ] **Step 2: Replace absolute spike totals with rate-normalized smoke assertions**

```python
def _check_ping(name, t_final_run=200.0, min_e_hz=5.0, min_i_hz=20.0):
    py = load_python_port(ROOT / "python" / PYTHON_BASE / name / "main.py")
    t_e, i_e, t_i, i_i, v_plot = py.run_smoke(t_final_run=t_final_run)
    e_hz = len(t_e) / py.num_e / t_final_run * 1000.0
    i_hz = len(t_i) / py.num_i / t_final_run * 1000.0
    assert e_hz > min_e_hz
    assert i_hz > min_i_hz
    assert v_plot.max() > 0.0
```

Call it for `M_CURRENT_PING_1`, `M_CURRENT_PING_1_CLOSEUP`, `M_CURRENT_PING_2_CLOSEUP`, and `M_CURRENT_PING_3_CLOSEUP`. For from-rest, run `simulate_from_rest(rest_duration=200.0, driven_duration=200.0)`, assert `t_e.size > 0`, `t_i.size > 0`, `t_e.min() > 1.0`, and rate-normalized firing above 5Hz E/20Hz I.

- [ ] **Step 3: Run the five exact tests and enforce the after-times**

Run:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_1_structure \
  tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_1_closeup_structure \
  tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_1_from_rest_structure \
  tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_2_closeup_structure \
  tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_3_closeup_structure \
  --durations=5
```

Expected: `5 passed`; every call below 30s; total below 115s; firing-rate and voltage assertions pass.

- [ ] **Step 4: Commit the chapter-32 smoke entry points**

```bash
git add \
  python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_1/main.py \
  python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_1_CLOSEUP/main.py \
  python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_1_FROM_REST/main.py \
  python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_2_CLOSEUP/main.py \
  python/32_M_Current_PING_and_Poisson_PING/M_CURRENT_PING_3_CLOSEUP/main.py \
  tests/test_ch32_m_current_ping_poisson_ping.py
git commit -m "test: bound chapter 32 m-current smoke runs"
```

---

### Task 6: JIT the Dense Chapter-35/36 f-I Scan Kernels

**Measured baselines and targets:** Both `test_rtm_f_i_curve_with_inhibition*` tests exceeded 115s. `test_rtm_f_i_curve_pulsed_excitation_structure` consumed more than 113s after two quick tests, and variant 2 exceeded 115s alone. These tests require 101-point or 201-point grids and, for pulsed excitation, more than ten samples on the 40Hz plateau; reducing to a smoke grid would remove that numerical contract. Required after-times: **< 30s per exact test**, including first-call compilation, and **< 90s for all four tests**.

**Files:**
- Modify: `python/35_Periodic_Inhibition/RTM_F_I_CURVE_WITH_INHIBITION/main.py`
- Modify: `python/35_Periodic_Inhibition/RTM_F_I_CURVE_WITH_INHIBITION_2/main.py`
- Modify: `python/36_F_I_Curves_Pulsed_Excitation/RTM_F_I_CURVE_PULSED_EXCITATION/main.py`
- Modify: `python/36_F_I_Curves_Pulsed_Excitation/RTM_F_I_CURVE_PULSED_EXCITATION_2/main.py`
- Test: `tests/test_ch35_periodic_inhibition.py:64-79`
- Test: `tests/test_ch36_f_i_curves_pulsed_excitation.py:17-38`

**Interfaces:**
- Preserve no-argument behavior and full-grid return values while adding optional current arrays: `f_i_curve_tonic(i_ext_values=i_ext_vec)`, `f_i_curve_periodic(g_amplitude=g_bar, i_ext_values=i_ext_vec)`, `f_i_curve_constant(i_ext_values=i_ext_vec)`, and `f_i_curve_pulsed(i_ext_values=i_ext_vec)`.
- Preserve the existing state carry-over between adjacent current values.

- [ ] **Step 1: Add failing cold-process performance guards**

Do not put wall-clock assertions inside pytest. Add the existing numerical assertions below before implementation and use `--durations` as the performance gate:

```python
assert len(py.f_vec_tonic) == 101
assert len(py.f_vec_periodic) == 101
assert py.f_vec_tonic[0] == 0.0
assert py.f_vec_tonic[-1] > 0.0
assert py.f_vec_periodic[-1] > 0.0

assert len(py.f_vec_constant) == 201
assert len(py.f_vec_pulsed) == 201
assert py.f_vec_constant[0] == 0.0
assert py.f_vec_constant[-1] > 0.0
assert py.f_vec_pulsed[-1] > 0.0
assert np.count_nonzero(np.round(py.f_vec_pulsed) == 40.0) > 10
```

- [ ] **Step 2: Parameterize the current grid and compile the dense loops**

In each of the four scan functions, add the `i_ext_values` parameter shown in Interfaces, replace `np.zeros(len(i_ext_vec))` with `np.zeros(len(i_ext_values))`, and replace `enumerate(i_ext_vec)` with `enumerate(i_ext_values)`. No equation or stopping condition changes.

At module scope after the current pure-Python definitions and before the first full-grid call, rebind the numerical call graph to cached dispatchers. This avoids copying the long, already-tested midpoint loop and makes the pure-Python outer-loop aliases available for three-point equivalence checks:

```python
from numba import njit

# Chapter 35 scripts
_f_i_curve_tonic_python = f_i_curve_tonic
_f_i_curve_periodic_python = f_i_curve_periodic
m_inf = njit(cache=True)(m_inf)
alpha_h = njit(cache=True)(alpha_h)
beta_h = njit(cache=True)(beta_h)
alpha_n = njit(cache=True)(alpha_n)
beta_n = njit(cache=True)(beta_n)
step = njit(cache=True)(step)
g_periodic = njit(cache=True)(g_periodic)
f_i_curve_tonic = njit(cache=True)(f_i_curve_tonic)
f_i_curve_periodic = njit(cache=True)(f_i_curve_periodic)

# Chapter 36 scripts
_f_i_curve_constant_python = f_i_curve_constant
_f_i_curve_pulsed_python = f_i_curve_pulsed
m_inf = njit(cache=True)(m_inf)
alpha_h = njit(cache=True)(alpha_h)
beta_h = njit(cache=True)(beta_h)
alpha_n = njit(cache=True)(alpha_n)
beta_n = njit(cache=True)(beta_n)
step = njit(cache=True)(step)
f_i_curve_constant = njit(cache=True)(f_i_curve_constant)
f_i_curve_pulsed = njit(cache=True)(f_i_curve_pulsed)
```

Do not JIT plotting code. Preserve current ordering because tonic scans carry the final state from one current to the next.

- [ ] **Step 3: Add direct kernel-versus-Python numerical equivalence tests**

For each family, retain a private pure-Python reference callable and compare a three-current slice:

```python
probe = np.array([i_ext_vec[0], i_ext_vec[len(i_ext_vec) // 2], i_ext_vec[-1]])
actual = f_i_curve_periodic(g_amplitude=g_bar, i_ext_values=probe)
expected = _f_i_curve_periodic_python(g_amplitude=g_bar, i_ext_values=probe)
np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-10)
```

Use the corresponding tonic/constant/pulsed function names and the same `rtol=0.0, atol=1e-10` threshold. Run these equivalence checks in the existing chapter-35/36 test files.

- [ ] **Step 4: Run the four exact tests twice**

Run twice to cover cold compilation and cache reuse:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_structure \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_2_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_2_structure \
  --durations=4
```

Expected on both runs: `4 passed`; each call below 30s; total below 90s; full 101/201-point arrays and plateau assertions unchanged.

- [ ] **Step 5: Commit the dense-scan kernels**

```bash
git add \
  python/35_Periodic_Inhibition/RTM_F_I_CURVE_WITH_INHIBITION/main.py \
  python/35_Periodic_Inhibition/RTM_F_I_CURVE_WITH_INHIBITION_2/main.py \
  python/36_F_I_Curves_Pulsed_Excitation/RTM_F_I_CURVE_PULSED_EXCITATION/main.py \
  python/36_F_I_Curves_Pulsed_Excitation/RTM_F_I_CURVE_PULSED_EXCITATION_2/main.py \
  tests/test_ch35_periodic_inhibition.py \
  tests/test_ch36_f_i_curves_pulsed_excitation.py
git commit -m "perf: jit dense non-brian f-i scans"
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
  tests/test_ch32_m_current_ping_poisson_ping.py \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_structure \
  tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_2_structure \
  tests/test_ch36_f_i_curves_pulsed_excitation.py \
  --durations=30
```

Expected: all selected tests pass; no call exceeds 45s; every task-specific after-time above is met.

- [ ] **Step 3: Run the remaining non-Brian suite in chapter batches capped at 115 seconds**

Use the same chapter/file batches recorded in `task-3-report.md`. Every command must use `timeout --signal=INT 115s`; if a batch reaches the cap, run its active exact node separately and record it before proceeding to later nodes. The known unrelated documentation inventory failure must be reported separately rather than hidden.

- [ ] **Step 4: Commit any test-only threshold corrections justified by the same numerical contracts**

```bash
git status --short
git diff --check
git log -7 --oneline
```

Expected: only files named in Tasks 1-6 are changed; no Brian file appears; `git diff --check` emits no output.
