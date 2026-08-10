# Chapter 20 Brian2 Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete and verify the four missing Chapter 20 Brian2 examples in the existing chapter notebook.

**Architecture:** Extend `brian/chapter20.ipynb` with callable calculation and simulation functions plus four plotting sections. Tests load notebook definitions with `load_notebook_as_module` and compare their observable outputs with the already MATLAB-validated Python ports, keeping MATLAB out of the normal validation loop.

**Tech Stack:** Python 3, Brian2, NumPy, Matplotlib, pytest, Jupyter notebook JSON

## Global Constraints

- Work on `develop_python`; do not create a separate worktree.
- Do not modify or stage existing changes in chapters 19, 33, or 34.
- Do not change verified Python or MATLAB implementations.
- Keep Chapter 20 in its existing single-notebook structure.
- Use test-driven development: observe each focused test fail before adding its implementation.
- An example is complete only after numeric and visual checks pass.
- Use the verified Python ports as the primary numeric oracle.

---

### Task 1: Jahr-Stevens Magnesium Block

**Files:**
- Modify: `brian/chapter20.ipynb`
- Create: `tests/test_ch20_brian_b_jahr_stevens.py`

**Interfaces:**
- Produces: `jahr_stevens_block(v_post) -> np.ndarray`, accepting scalar or array-like membrane voltages in mV.

- [ ] **Step 1: Write the failing behavior test**

The production mutation caught by this test is an incorrect voltage coefficient, magnesium divisor, or sign in the block equation.

```python
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]


def test_brian_b_jahr_stevens_matches_verified_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter20.ipynb")
    py = load_python_port(
        ROOT / "python" / "20_Chemical_Synapses" / "B_JAHR_STEVENS" / "main.py"
    )
    np.testing.assert_allclose(ns.jahr_stevens_block(py.v), py.B, atol=1e-12)
```

- [ ] **Step 2: Run the focused test and verify RED**

Run: `python3 -m pytest tests/test_ch20_brian_b_jahr_stevens.py -v`

Expected: FAIL because `jahr_stevens_block` is absent.

- [ ] **Step 3: Add the minimal notebook section**

Add a `### 20.x Voltage-Dependent Magnesium Block` markdown cell, followed by a definition cell containing:

```python
def jahr_stevens_block(v_post):
    v_post = np.asarray(v_post)
    return 1.0 / (1.0 + np.exp(-0.062 * v_post) / 3.57)
```

Add a plotting cell using `v_post = np.arange(-100, 51)` with the same labels and limits as the Python figure.

- [ ] **Step 4: Run the focused test and verify GREEN**

Run: `python3 -m pytest tests/test_ch20_brian_b_jahr_stevens.py -v`

Expected: PASS.

- [ ] **Step 5: Commit only Task 1 files**

```bash
git add brian/chapter20.ipynb tests/test_ch20_brian_b_jahr_stevens.py
git commit -m "feat: add Chapter 20 Jahr-Stevens Brian example"
```

### Task 2: RTM Plot Q

**Files:**
- Modify: `brian/chapter20.ipynb`
- Create: `tests/test_ch20_brian_rtm_plot_q.py`

**Interfaces:**
- Consumes: existing `simulate_RTM_neuron_q_s(...)` from `brian/chapter20.ipynb`.
- Produces: the same function with `StateMonitor` variables `vm`, `q`, and `s`.

- [ ] **Step 1: Write the failing behavior test**

The production mutation caught by this test is omitting `q` from monitoring or using the wrong release/rise equation.

```python
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_brian_rtm_plot_q_matches_verified_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter20.ipynb")
    py = load_python_port(
        ROOT / "python" / "20_Chemical_Synapses" / "RTM_PLOT_Q" / "main.py"
    )
    py.tau_r, py.tau_d = 0.1, 2.0
    py.t = np.arange(0.0, 100.0, 0.01)
    py_trace = py.odeint(py.derivative, py.initial_condition(-70.0), py.t)
    sm = ns.simulate_RTM_neuron_q_s(
        1.0 * ns.b2.uA,
        tau_r=0.1 * ns.b2.ms,
        tau_d=2.0 * ns.b2.ms,
        tau_d_q=2.0 * ns.b2.ms,
        simulation_time=100 * ns.b2.ms,
    )
    assert trace_rmse(py.t, py_trace[:, 0], sm.t / ns.b2.ms, sm.vm[0] / ns.b2.mV) < 2.0
    assert trace_rmse(py.t, py_trace[:, 3], sm.t / ns.b2.ms, sm.q[0]) < 0.03
    assert trace_rmse(py.t, py_trace[:, 4], sm.t / ns.b2.ms, sm.s[0]) < 0.03
```

- [ ] **Step 2: Run the focused test and verify RED**

Run: `python3 -m pytest tests/test_ch20_brian_rtm_plot_q.py -v`

Expected: FAIL because the returned monitor has no `q` state.

- [ ] **Step 3: Make the minimal monitor and plotting changes**

Change the state monitor inside `simulate_RTM_neuron_q_s` to:

```python
st_mon = b2.StateMonitor(neurons, ["vm", "q", "s"], record=True)
```

Add an `RTM_PLOT_Q` plotting cell that calls the function with the parameters in the test and plots voltage in the first panel and `s`/`q` in the second panel.

- [ ] **Step 4: Run focused and existing Chapter 20 tests**

Run: `python3 -m pytest tests/test_ch20_brian_rtm_plot_q.py tests/test_ch20_brian_rtm_plot_s_two_variables.py tests/test_ch20_brian_s_buildup.py tests/test_ch20_brian_s_slow_buildup.py -v`

Expected: all PASS.

- [ ] **Step 5: Commit only Task 2 files**

```bash
git add brian/chapter20.ipynb tests/test_ch20_brian_rtm_plot_q.py
git commit -m "feat: add Chapter 20 RTM q-gate Brian example"
```

### Task 3: Prescribed Synaptic Peak Timing

**Files:**
- Modify: `brian/chapter20.ipynb`
- Create: `tests/test_ch20_brian_rtm_plot_s_prescribe_tau_peak.py`

**Interfaces:**
- Produces: `synaptic_peak_time(tau_d: float, tau_r: float, tau_d_q: float, dt: float = 0.01) -> float` in ms.
- Produces: `release_time_constant(tau_d: float, tau_r: float, tau_peak: float) -> float` in ms.
- Consumes: `simulate_RTM_neuron_q_s(...)`.

- [ ] **Step 1: Write the failing helper and trace test**

The production mutations caught are an incorrect decay term in the peak integrator, reversed bisection branch, or failure to pass the solved `tau_d_q` into Brian2.

```python
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_brian_prescribed_tau_peak_matches_verified_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter20.ipynb")
    py = load_python_port(
        ROOT / "python" / "20_Chemical_Synapses"
        / "RTM_PLOT_S_PRESCRIBE_TAU_PEAK" / "main.py"
    )
    tau_d_q = ns.release_time_constant(300.0, 10.0, 20.0)
    expected = py.tau_d_q_function(300.0, 10.0, 20.0)
    assert np.isclose(tau_d_q, expected, rtol=1e-8)
    assert np.isclose(ns.synaptic_peak_time(300.0, 10.0, tau_d_q), 20.0, atol=0.02)

    py.tau_d, py.tau_r, py.tau_d_q = 300.0, 10.0, expected
    py.i_ext = 0.12
    py.t = np.arange(0.0, 500.0, 0.01)
    py_trace = py.odeint(py.derivative, py.initial_condition(-70.0), py.t)
    sm = ns.simulate_RTM_neuron_q_s(
        0.12 * ns.b2.uA,
        tau_r=10 * ns.b2.ms,
        tau_d=300 * ns.b2.ms,
        tau_d_q=tau_d_q * ns.b2.ms,
        simulation_time=500 * ns.b2.ms,
    )
    assert trace_rmse(py.t, py_trace[:, 4], sm.t / ns.b2.ms, sm.s[0]) < 0.03
```

- [ ] **Step 2: Run the focused test and verify RED**

Run: `python3 -m pytest tests/test_ch20_brian_rtm_plot_s_prescribe_tau_peak.py -v`

Expected: FAIL because both peak-timing helpers are absent.

- [ ] **Step 3: Implement the peak-time and bisection helpers**

Translate the Python reference's midpoint peak search exactly, including `-s / tau_d` in every derivative evaluation. Implement bisection with a `1e-12` ms interval tolerance and return plain float milliseconds.

- [ ] **Step 4: Add the plotting cell**

Calculate `(tau_r, tau_peak) = (10, 20)` and `(100, 150)` with `tau_d = 300` ms, run both 2000 ms simulations, and plot voltage plus the two `s` traces in three stacked panels.

- [ ] **Step 5: Run focused and dependent tests**

Run: `python3 -m pytest tests/test_ch20_brian_rtm_plot_s_prescribe_tau_peak.py tests/test_ch20_brian_rtm_plot_s_two_variables.py -v`

Expected: all PASS.

- [ ] **Step 6: Commit only Task 3 files**

```bash
git add brian/chapter20.ipynb tests/test_ch20_brian_rtm_plot_s_prescribe_tau_peak.py
git commit -m "feat: prescribe Chapter 20 Brian synaptic peak timing"
```

### Task 4: Autaptic RTM F-I Curve

**Files:**
- Modify: `brian/chapter20.ipynb`
- Create: `tests/test_ch20_brian_rtm_with_autapse_f_i_curve.py`

**Interfaces:**
- Produces: `simulate_rtm_autapse_frequency(i_ext, initial_state=None, dt=0.005*b2.ms) -> tuple[float, dict[str, float]]`.
- Produces: `sweep_rtm_autapse(i_ext_values=np.linspace(0, 0.15, 31)) -> dict` with keys `i_ext`, `f_forward`, `f_backward`, `i_c`, and `i_star`.
- Consumes: `release_time_constant(5.0, 0.2, 0.6)` from Task 3.

- [ ] **Step 1: Write the failing observable-behavior test**

The production mutations caught are a missing autaptic conductance, wrong synaptic reversal potential, loss of state continuation between sweep points, or a reversed backward sweep.

```python
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]


def test_brian_rtm_autapse_f_i_matches_verified_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter20.ipynb")
    py = load_python_port(
        ROOT / "python" / "20_Chemical_Synapses"
        / "RTM_WITH_AUTAPSE_F_I_CURVE" / "main.py"
    )
    result = ns.sweep_rtm_autapse(py.i_ext_vec)
    np.testing.assert_allclose(result["f_backward"], py.f_backward, atol=1.0)
    assert np.isclose(result["i_c"], py.I_c, atol=0.0026)
    assert np.isclose(result["i_star"], py.I_star, atol=0.0026)
```

- [ ] **Step 2: Run the focused test and verify RED**

Run: `python3 -m pytest tests/test_ch20_brian_rtm_with_autapse_f_i_curve.py -v`

Expected: FAIL because `sweep_rtm_autapse` is absent.

- [ ] **Step 3: Implement one-current Brian2 simulation**

Use the existing RTM voltage/gating equations, add `g_syn * s * (0*mV - vm)` with `g_syn = 0.1*mS`, and use the `q/s` kinetics with `tau_d=5 ms`, `tau_r=0.2 ms`, and solved `tau_d_q`. Detect downward `-20 mV` crossings from a state monitor, return zero after the same steady-state criterion as the Python port, and otherwise calculate Hz from the third-to-fourth spike interval. Return the final `vm`, `h`, `n`, `q`, and `s` so the next current starts from the previous state.

- [ ] **Step 4: Implement the forward/backward sweep**

Iterate currents in ascending order, retain state, then iterate the same values in descending order without resetting state. Compute `i_c` and `i_star` as midpoints immediately above the last zero-frequency point in the forward and backward arrays.

- [ ] **Step 5: Add the full plotting cell**

Run 31 currents from `0.0` through `0.15` uA, plot filled forward and open backward markers, and draw the two critical-current lines matching the verified Python figure.

- [ ] **Step 6: Run the focused test and verify GREEN**

Run: `python3 -m pytest tests/test_ch20_brian_rtm_with_autapse_f_i_curve.py -v`

Expected: PASS. If runtime exceeds five minutes, profile the real simulation and reduce repeated monitor allocation without replacing it with mocked behavior.

- [ ] **Step 7: Commit only Task 4 files**

```bash
git add brian/chapter20.ipynb tests/test_ch20_brian_rtm_with_autapse_f_i_curve.py
git commit -m "feat: add Chapter 20 Brian autapse F-I curve"
```

### Task 5: Figures, Progress, and Full Verification

**Files:**
- Modify: `brian/chapter20.ipynb`
- Modify: `PROGRESS.md`

**Interfaces:**
- Consumes: all Task 1-4 notebook functions and plotting cells.
- Produces: four visually checked notebook figures and Chapter 20 Brian status `8/8`.

- [ ] **Step 1: Execute the Chapter 20 notebook headlessly**

Run: `MPLBACKEND=Agg jupyter nbconvert --to notebook --execute brian/chapter20.ipynb --output /tmp/chapter20-executed.ipynb --ExecutePreprocessor.timeout=900`

Expected: execution completes without cell errors. Keep the tracked notebook output-free according to the repository's `nbstripout` policy.

- [ ] **Step 2: Visually inspect the four new figures**

Compare the executed notebook outputs with:

- `python/20_Chemical_Synapses/B_JAHR_STEVENS/fig.png`
- `python/20_Chemical_Synapses/RTM_PLOT_Q/fig_20_3.png`
- `python/20_Chemical_Synapses/RTM_PLOT_S_PRESCRIBE_TAU_PEAK/fig_20_5.png`
- `python/20_Chemical_Synapses/RTM_WITH_AUTAPSE_F_I_CURVE/fig.png`

Check curve direction, panel count, axes, trace timing, forward/backward markers, and critical-current positions.

- [ ] **Step 3: Update progress only after visual verification**

In `PROGRESS.md`, change the Brian cells for the four named Chapter 20 examples from `[ ]` to `[x]`, and change the total from `brian 46/256` to `brian 50/256`. Leave Python and n/a totals unchanged.

- [ ] **Step 4: Run all Chapter 20 tests**

Run: `MPLBACKEND=Agg python3 -m pytest tests/test_ch20*.py -v`

Expected: all PASS.

- [ ] **Step 5: Run the complete suite**

Run: `MPLBACKEND=Agg python3 -m pytest -q`

Expected: all PASS.

- [ ] **Step 6: Audit scope**

Run: `git status --short`

Run: `git diff -- PROGRESS.md brian/chapter20.ipynb tests/test_ch20_brian_b_jahr_stevens.py tests/test_ch20_brian_rtm_plot_q.py tests/test_ch20_brian_rtm_plot_s_prescribe_tau_peak.py tests/test_ch20_brian_rtm_with_autapse_f_i_curve.py`

Confirm the diff contains no chapter 19, 33, or 34 changes.

- [ ] **Step 7: Commit final progress and any output-stripping normalization**

```bash
git add PROGRESS.md brian/chapter20.ipynb
git commit -m "docs: mark Chapter 20 Brian2 complete"
```
