# Chapter 31 Brian2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement and verify the 14 Chapter 31 examples that genuinely use neuronal simulation in Brian2, while marking the two analytical phase-map examples not applicable.

**Architecture:** Create one `brian/chapter31.ipynb` with reusable WB and RTM equations, deterministic input/connectivity builders, one parameterized ING-network simulator, and one E/I entrainment simulator. Keep plotting behind `if __name__ == "__main__"` guards so tests can load the complete notebook without running full examples.

**Tech Stack:** Python 3, Brian2, NumPy, Matplotlib, pytest, Jupyter notebook JSON

## Global Constraints

- Preserve all existing uncommitted work; do not modify or stage unrelated Chapter 19, 20, 34, or 35 files.
- Do not modify the verified files under `python/31_ING_Rhythms` or `matlab/31`.
- Implement only `1_CELL_ING`, `1_CELL_ING_CONDITION_NUMBERS`, `ING_1` through `ING_10`, `ING_ENTRAINING_E_CELLS`, and `ING_ENTRAINING_E_CELLS_2`.
- Do not copy `ABSTRACT_PULSE_COUPLING_INH` or `ABSTRACT_PULSE_COUPLING_INH_2` into the Brian notebook; mark both `n/a` in `PROGRESS.md`.
- Preserve the reference timestep of `0.01 ms` by default and make duration and population sizes overrideable for focused tests.
- Generate stochastic drives, connectivity, and phase samples with NumPy and assign the resulting arrays to Brian2 explicitly.
- Call `brian2.start_scope()` at the start of every public simulation.
- Use test-driven development: observe each focused test fail before adding the corresponding notebook definitions.
- A simulation example is complete only after numeric tests and full-duration visual checks pass.

---

### Task 1: Shared Notebook Foundation and WB Kinetics

**Files:**
- Create: `brian/chapter31.ipynb`
- Create: `tests/test_ch31_brian_shared.py`

**Interfaces:**
- Produces: `tau_peak(tau_d_ms, tau_r_ms, tau_dq_ms, dt_ms=0.01) -> float`.
- Produces: `solve_tau_dq(tau_d_ms, tau_r_ms, tau_peak_ms) -> float`.
- Produces: `validate_probability(name, value) -> float`.
- Produces: `WB_EQS`, a Brian2 equation string used by later tasks.
- Produces: `initial_wb_state(count, mode, rng) -> dict[str, np.ndarray]`.

- [ ] **Step 1: Write the failing shared-helper tests**

```python
from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "brian" / "chapter31.ipynb"


def test_ch31_notebook_exposes_shared_helpers():
    ns = load_notebook_as_module(NOTEBOOK)
    py = load_python_port(
        ROOT / "python" / "31_ING_Rhythms" / "1_CELL_ING_CONDITION_NUMBERS" / "main.py"
    )
    expected = py.tau_d_q_function(9.0, 0.5, 0.5)
    actual = ns.solve_tau_dq(9.0, 0.5, 0.5)
    assert np.isclose(actual, expected, rtol=1e-10)
    assert np.isclose(ns.tau_peak(9.0, 0.5, actual), 0.5, atol=0.02)


def test_probability_validation_is_explicit():
    ns = load_notebook_as_module(NOTEBOOK)
    assert ns.validate_probability("p_ii", 0.5) == 0.5
    with pytest.raises(ValueError, match="p_ii"):
        ns.validate_probability("p_ii", -0.01)
    with pytest.raises(ValueError, match="p_gap"):
        ns.validate_probability("p_gap", 1.01)
```

- [ ] **Step 2: Run the test and verify RED**

Run: `python3 -m pytest tests/test_ch31_brian_shared.py -v`

Expected: FAIL because `brian/chapter31.ipynb` does not exist.

- [ ] **Step 3: Create the notebook title/import cells and numerical helpers**

The first code cell must contain:

```python
import brian2 as b2
import matplotlib.pyplot as plt
import numpy as np


def validate_probability(name, value):
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must be in [0, 1], got {value}")
    return value


def tau_peak(tau_d_ms, tau_r_ms, tau_dq_ms, dt_ms=0.01):
    s = 0.0
    t = 0.0
    ds = np.exp(-t / tau_dq_ms) * (1.0 - s) / tau_r_ms - s / tau_d_ms
    while ds > 0.0:
        t_old, ds_old = t, ds
        s_mid = s + 0.5 * dt_ms * ds
        ds_mid = (
            np.exp(-(t + 0.5 * dt_ms) / tau_dq_ms)
            * (1.0 - s_mid) / tau_r_ms
            - s_mid / tau_d_ms
        )
        s += dt_ms * ds_mid
        t += dt_ms
        ds = np.exp(-t / tau_dq_ms) * (1.0 - s) / tau_r_ms - s / tau_d_ms
    return (t_old * (-ds) + t * ds_old) / (ds_old - ds)


def solve_tau_dq(tau_d_ms, tau_r_ms, tau_peak_ms):
    left = 1.0
    while tau_peak(tau_d_ms, tau_r_ms, left) > tau_peak_ms:
        left *= 0.5
    right = tau_r_ms
    while tau_peak(tau_d_ms, tau_r_ms, right) < tau_peak_ms:
        right *= 2.0
    while right - left > 1e-12:
        middle = 0.5 * (left + right)
        if tau_peak(tau_d_ms, tau_r_ms, middle) <= tau_peak_ms:
            left = middle
        else:
            right = middle
    return 0.5 * (left + right)
```

- [ ] **Step 4: Add the shared WB Brian2 equation cell**

Use volt, second, amp, siemens, and farad consistently. Define the following equation string exactly once:

```python
WB_EQS = """
dv/dt = (g_l*(E_l-v) + g_k*n**4*(E_k-v)
         + g_na*m_inf**3*h*(E_na-v) + i_ext + I_chem + I_gap)/C : volt
dh/dt = (h_inf-h)/tau_h : 1
dn/dt = (n_inf-n)/tau_n : 1
dq/dt = 0.5*(1+tanh(v/(10*mV)))*(1-q)/(0.1*ms) - q/tau_dq : 1
ds/dt = q*(1-s)/tau_r - s/tau_d : 1
m_inf = alpha_m/(alpha_m+beta_m) : 1
h_inf = alpha_h/(alpha_h+beta_h) : 1
n_inf = alpha_n/(alpha_n+beta_n) : 1
alpha_m = 0.1/mV*(v+35*mV)/(1-exp(-(v+35*mV)/(10*mV)))/ms : Hz
beta_m = 4*exp(-(v+60*mV)/(18*mV))/ms : Hz
alpha_h = 0.07*exp(-(v+58*mV)/(20*mV))/ms : Hz
beta_h = 1/(exp(-(v+28*mV)/(10*mV))+1)/ms : Hz
alpha_n = -0.01/mV*(v+34*mV)/(exp(-(v+34*mV)/(10*mV))-1)/ms : Hz
beta_n = 0.125*exp(-(v+44*mV)/(80*mV))/ms : Hz
tau_h = 1/(5*(alpha_h+beta_h)) : second
tau_n = 1/(5*(alpha_n+beta_n)) : second
i_ext : amp
I_chem : amp
I_gap : amp
tau_dq : second
tau_r : second
tau_d : second
C : farad (constant)
g_l : siemens (constant)
g_k : siemens (constant)
g_na : siemens (constant)
E_l : volt (constant)
E_k : volt (constant)
E_na : volt (constant)
"""
```

Define `initial_wb_state` with `uniform_random` producing `v` uniformly in
`[-100, 50] mV` and `h/n` uniformly in `[0, 1]`; define `fixed` as
`v=-75 mV`, `h=n=0.1`; both modes initialize `q=s=0`. Reject other mode names
with `ValueError`.

- [ ] **Step 5: Run the focused tests and verify GREEN**

Run: `python3 -m pytest tests/test_ch31_brian_shared.py -v`

Expected: 2 passed.

- [ ] **Step 6: Commit Task 1 only**

```bash
git add brian/chapter31.ipynb tests/test_ch31_brian_shared.py
git commit -m "feat: add Chapter 31 Brian shared WB model"
```

---

### Task 2: Single-Cell ING Simulation

**Files:**
- Modify: `brian/chapter31.ipynb`
- Create: `tests/test_ch31_brian_1_cell_ing.py`

**Interfaces:**
- Consumes: `WB_EQS`, `solve_tau_dq(...)`.
- Produces: `simulate_single_cell_ing(i_ext=1.5*b2.uA, g_ii=0.5*b2.msiemens, tau_d=9*b2.ms, duration=200*b2.ms, dt=0.01*b2.ms, record=True) -> dict`.
- Result keys: `t_ms`, `v_mv`, `spike_ms`, `period_ms`, and `frequency_hz`.

- [ ] **Step 1: Write the failing behavior test**

```python
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port, spike_times, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_brian_single_cell_ing_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    py = load_python_port(ROOT / "python" / "31_ING_Rhythms" / "1_CELL_ING" / "main.py")
    py.tau_dq_i = py.tau_d_q_function(py.tau_d_i, py.tau_r_i, py.tau_peak_i)
    t = np.arange(0.0, py.t_final, py.dt)
    sol = py.odeint(py.derivative, py.x0, t)
    expected_spikes = spike_times(t, sol[:, 0], threshold=-20.0)

    result = ns.simulate_single_cell_ing()
    actual_spikes = np.asarray(result["spike_ms"])

    assert trace_rmse(t, sol[:, 0], result["t_ms"], result["v_mv"]) < 3.0
    assert len(actual_spikes) == len(expected_spikes)
    assert np.isclose(result["period_ms"], np.diff(expected_spikes)[-1], atol=0.5)
```

- [ ] **Step 2: Run the test and verify RED**

Run: `python3 -m pytest tests/test_ch31_brian_1_cell_ing.py -v`

Expected: FAIL because `simulate_single_cell_ing` is absent.

- [ ] **Step 3: Implement the one-cell Brian simulation**

Create one `NeuronGroup` with `WB_EQS`, method `rk2`, and threshold
`v > -20*mV`. Assign the reference constants:

```python
C = 1.0 * b2.ufarad
g_l = 0.1 * b2.msiemens
g_k = 9.0 * b2.msiemens
g_na = 35.0 * b2.msiemens
E_l = -65.0 * b2.mV
E_k = -90.0 * b2.mV
E_na = 55.0 * b2.mV
```

Implement the autaptic current with a one-edge `Synapses` object:

```python
chem = b2.Synapses(
    cells,
    cells,
    model="g : siemens\nI_chem_post = g*s_pre*(-75*mV-v_post) : amp (summed)",
)
chem.connect(i=[0], j=[0])
chem.g = g_ii
```

Initialize `v=-75 mV`, `h=n=0.1`, `q=s=0`; calculate `tau_dq` with
`solve_tau_dq`; record `v`; run; and return plain NumPy arrays in milliseconds
and millivolts. Calculate period from the last two spike-monitor times and
raise `RuntimeError("fewer than two spikes")` if a parameter set does not
produce two spikes.

- [ ] **Step 4: Run the focused and shared tests**

Run: `python3 -m pytest tests/test_ch31_brian_1_cell_ing.py tests/test_ch31_brian_shared.py -v`

Expected: 3 passed.

- [ ] **Step 5: Add the guarded reference plot**

Add a markdown heading for `1_CELL_ING`, then:

```python
if __name__ == "__main__":
    one_cell = simulate_single_cell_ing()
    fig, ax = plt.subplots(figsize=(7, 3))
    ax.plot(one_cell["t_ms"], one_cell["v_mv"], color="k", lw=2)
    ax.set(xlabel="time [ms]", ylabel="v [mV]", ylim=(-100, 50))
    fig.tight_layout()
```

- [ ] **Step 6: Commit Task 2 only**

```bash
git add brian/chapter31.ipynb tests/test_ch31_brian_1_cell_ing.py
git commit -m "feat: add Chapter 31 single-cell ING Brian simulation"
```

---

### Task 3: Single-Cell ING Condition Numbers

**Files:**
- Modify: `brian/chapter31.ipynb`
- Create: `tests/test_ch31_brian_1_cell_ing_condition_numbers.py`

**Interfaces:**
- Consumes: `simulate_single_cell_ing(...)`.
- Produces: `compute_ing_condition_numbers(duration=200*b2.ms) -> dict`.
- Result keys: `base_period_ms`, `pct_i_ext`, `pct_g_ii`, `pct_tau_d`.

- [ ] **Step 1: Write the failing test**

```python
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]


def test_brian_ing_condition_numbers_match_verified_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    py = load_python_port(
        ROOT / "python" / "31_ING_Rhythms" / "1_CELL_ING_CONDITION_NUMBERS" / "main.py"
    )
    result = ns.compute_ing_condition_numbers()
    assert np.isclose(result["pct_i_ext"], py.pct_i_ext, atol=0.08)
    assert np.isclose(result["pct_g_ii"], py.pct_g_ii, atol=0.08)
    assert np.isclose(result["pct_tau_d"], py.pct_tau_d, atol=0.08)
```

- [ ] **Step 2: Run the test and verify RED**

Run: `python3 -m pytest tests/test_ch31_brian_1_cell_ing_condition_numbers.py -v`

Expected: FAIL because `compute_ing_condition_numbers` is absent.

- [ ] **Step 3: Implement the four-run sensitivity calculation**

```python
def compute_ing_condition_numbers(duration=200*b2.ms):
    base = simulate_single_cell_ing(duration=duration, record=False)
    reduced_i = simulate_single_cell_ing(i_ext=1.5*0.99*b2.uA, duration=duration, record=False)
    raised_g = simulate_single_cell_ing(g_ii=0.5*1.01*b2.msiemens, duration=duration, record=False)
    raised_tau = simulate_single_cell_ing(tau_d=9*1.01*b2.ms, duration=duration, record=False)
    period = base["period_ms"]
    return {
        "base_period_ms": period,
        "pct_i_ext": 100.0*(period-reduced_i["period_ms"])/period,
        "pct_g_ii": 100.0*(period-raised_g["period_ms"])/period,
        "pct_tau_d": 100.0*(period-raised_tau["period_ms"])/period,
    }
```

Add a guarded cell that prints the four returned values using the same labels
as the Python port.

- [ ] **Step 4: Run all single-cell tests**

Run: `python3 -m pytest tests/test_ch31_brian_shared.py tests/test_ch31_brian_1_cell_ing.py tests/test_ch31_brian_1_cell_ing_condition_numbers.py -v`

Expected: 4 passed.

- [ ] **Step 5: Commit Task 3 only**

```bash
git add brian/chapter31.ipynb tests/test_ch31_brian_1_cell_ing_condition_numbers.py
git commit -m "feat: add Chapter 31 ING condition numbers"
```

---

### Task 4: ING Configuration and Deterministic Connectivity

**Files:**
- Modify: `brian/chapter31.ipynb`
- Create: `tests/test_ch31_brian_ing_configuration.py`

**Interfaces:**
- Produces: `ING_CONFIGS: dict[str, dict]`.
- Produces: `wb_phase_initial_state(i_ext_ua, phase, dt_ms=0.01) -> dict[str, np.ndarray]`.
- Produces: `build_ing_inputs(config, seed=None, num_i=None) -> dict`.
- Input result keys: `i_ext_ua`, `g_ii_ms`, `g_gap_ms`, `initial_state`.

- [ ] **Step 1: Write failing configuration tests**

```python
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_all_ten_ing_configurations_are_exact():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    expected = {
        "ING_1":  (0.00, 0.5, 1.00, 0.00, 1.00, False, "uniform_random"),
        "ING_2":  (0.03, 0.5, 1.00, 0.00, 1.00, False, "phase"),
        "ING_3":  (0.00, 0.5, 0.85, 0.00, 1.00, False, "phase"),
        "ING_4":  (0.00, 0.5, 0.85, 0.00, 1.00, True,  "phase"),
        "ING_5":  (0.05, 0.5, 0.50, 0.00, 1.00, False, "phase"),
        "ING_6":  (0.05, 0.5, 0.50, 0.10, 0.05, False, "phase"),
        "ING_7":  (0.00, 0.5, 1.00, 0.00, 1.00, False, "phase"),
        "ING_8":  (0.00, 0.5, 1.00, 0.00, 0.05, False, "phase"),
        "ING_9":  (0.05, 0.5, 1.00, 0.00, 0.05, False, "phase"),
        "ING_10": (0.05, 0.5, 1.00, 0.04, 0.05, False, "phase"),
    }
    actual = {
        name: (
            cfg["sigma_i"], cfg["g_hat_ii"], cfg["p_ii"],
            cfg["g_hat_gap"], cfg["p_gap"], cfg["fixed_indegree"], cfg["init_mode"],
        )
        for name, cfg in ns.ING_CONFIGS.items()
    }
    assert actual == expected


def test_connectivity_normalization_and_symmetry():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    built = ns.build_ing_inputs(ns.ING_CONFIGS["ING_6"], seed=63806, num_i=40)
    assert built["g_ii_ms"].shape == (40, 40)
    assert built["g_gap_ms"].shape == (40, 40)
    assert np.allclose(built["g_gap_ms"], built["g_gap_ms"].T)
    assert np.allclose(np.diag(built["g_gap_ms"]), 0.0)


def test_fixed_indegree_has_equal_column_counts():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    built = ns.build_ing_inputs(ns.ING_CONFIGS["ING_4"], seed=63806, num_i=40)
    counts = np.count_nonzero(built["g_ii_ms"], axis=0)
    assert np.array_equal(counts, np.full(40, round(0.85*40)))
```

- [ ] **Step 2: Run the test and verify RED**

Run: `python3 -m pytest tests/test_ch31_brian_ing_configuration.py -v`

Expected: FAIL because `ING_CONFIGS` and `build_ing_inputs` are absent.

- [ ] **Step 3: Add the exact configuration mapping**

```python
ING_CONFIGS = {
    "ING_1":  dict(sigma_i=0.00, g_hat_ii=0.5, p_ii=1.00, g_hat_gap=0.00, p_gap=1.00, fixed_indegree=False, init_mode="uniform_random", seed=124875),
    "ING_2":  dict(sigma_i=0.03, g_hat_ii=0.5, p_ii=1.00, g_hat_gap=0.00, p_gap=1.00, fixed_indegree=False, init_mode="phase", seed=63806),
    "ING_3":  dict(sigma_i=0.00, g_hat_ii=0.5, p_ii=0.85, g_hat_gap=0.00, p_gap=1.00, fixed_indegree=False, init_mode="phase", seed=63806),
    "ING_4":  dict(sigma_i=0.00, g_hat_ii=0.5, p_ii=0.85, g_hat_gap=0.00, p_gap=1.00, fixed_indegree=True,  init_mode="phase", seed=63806),
    "ING_5":  dict(sigma_i=0.05, g_hat_ii=0.5, p_ii=0.50, g_hat_gap=0.00, p_gap=1.00, fixed_indegree=False, init_mode="phase", seed=63806),
    "ING_6":  dict(sigma_i=0.05, g_hat_ii=0.5, p_ii=0.50, g_hat_gap=0.10, p_gap=0.05, fixed_indegree=False, init_mode="phase", seed=63806),
    "ING_7":  dict(sigma_i=0.00, g_hat_ii=0.5, p_ii=1.00, g_hat_gap=0.00, p_gap=1.00, fixed_indegree=False, init_mode="phase", seed=63806),
    "ING_8":  dict(sigma_i=0.00, g_hat_ii=0.5, p_ii=1.00, g_hat_gap=0.00, p_gap=0.05, fixed_indegree=False, init_mode="phase", seed=63806),
    "ING_9":  dict(sigma_i=0.05, g_hat_ii=0.5, p_ii=1.00, g_hat_gap=0.00, p_gap=0.05, fixed_indegree=False, init_mode="phase", seed=63806),
    "ING_10": dict(sigma_i=0.05, g_hat_ii=0.5, p_ii=1.00, g_hat_gap=0.04, p_gap=0.05, fixed_indegree=False, init_mode="phase", seed=63806),
}
```

- [ ] **Step 4: Implement deterministic inputs**

Add the following matrix construction; call `wb_phase_initial_state` in the
`phase` branch after generating `phase = rng.random(num_i)`:

```python
def build_ing_inputs(config, seed=None, num_i=None):
    config = dict(config)
    num_i = 100 if num_i is None else int(num_i)
    if num_i < 1:
        raise ValueError("num_i must be positive")
    p_ii = validate_probability("p_ii", config["p_ii"])
    p_gap = validate_probability("p_gap", config["p_gap"])
    if p_ii == 0.0 and config["g_hat_ii"] != 0.0:
        raise ValueError("nonzero g_hat_ii requires p_ii > 0")
    if p_gap == 0.0 and config["g_hat_gap"] != 0.0:
        raise ValueError("nonzero g_hat_gap requires p_gap > 0")

    rng = np.random.default_rng(config["seed"] if seed is None else seed)
    i_ext_ua = 1.5 * (1.0 + config["sigma_i"] * rng.standard_normal(num_i))

    g_ii_ms = np.zeros((num_i, num_i))
    if config["g_hat_ii"] != 0.0:
        weight = config["g_hat_ii"] / (num_i * p_ii)
        if config["fixed_indegree"]:
            g_ii_ms.fill(weight)
            omit = round(num_i - p_ii * num_i)
            for post in range(num_i):
                drop = rng.choice(num_i, size=omit, replace=False)
                g_ii_ms[drop, post] = 0.0
        else:
            g_ii_ms = weight * (rng.random((num_i, num_i)) < p_ii)

    g_gap_ms = np.zeros((num_i, num_i))
    if config["g_hat_gap"] != 0.0 and num_i > 1:
        weight = config["g_hat_gap"] / (p_gap * (num_i - 1))
        for left in range(num_i - 1):
            for right in range(left + 1, num_i):
                if rng.random() < p_gap:
                    g_gap_ms[left, right] = weight
                    g_gap_ms[right, left] = weight

    if config["init_mode"] == "phase":
        initial_state = wb_phase_initial_state(i_ext_ua, rng.random(num_i))
    else:
        initial_state = initial_wb_state(num_i, config["init_mode"], rng)
    return dict(
        i_ext_ua=i_ext_ua,
        g_ii_ms=g_ii_ms,
        g_gap_ms=g_gap_ms,
        initial_state=initial_state,
    )
```

`wb_phase_initial_state` must:

1. validate `p_ii` and `p_gap`;
2. reject nonzero strength with zero probability;
3. create `i_ext_ua = 1.5*(1 + sigma_i*rng.standard_normal(num_i))`;
4. create `g_ii_ms` with weight `g_hat_ii/(num_i*p_ii)`;
5. for fixed indegree, drop exactly `round(num_i-p_ii*num_i)` entries from each postsynaptic column without replacement;
6. create a symmetric, zero-diagonal `g_gap_ms` with pair weight `g_hat_gap/(p_gap*(num_i-1))`; and
7. return arrays for `v`, `h`, `n`, `q`, and `s`.

Implement `wb_phase_initial_state(i_ext_ua, phase, dt_ms=0.01)` by integrating
the uncoupled WB equations with the reference midpoint method, stopping after
the third downward `-20 mV` crossing of every cell, and interpolating each
cell between its second and third crossings at `second + phase*(third-second)`.
Return `q=s=0` with the interpolated `v/h/n` values.

- [ ] **Step 5: Run configuration and shared tests**

Run: `python3 -m pytest tests/test_ch31_brian_ing_configuration.py tests/test_ch31_brian_shared.py -v`

Expected: 5 passed.

- [ ] **Step 6: Commit Task 4 only**

```bash
git add brian/chapter31.ipynb tests/test_ch31_brian_ing_configuration.py
git commit -m "feat: define Chapter 31 ING configurations"
```

---

### Task 5: Parameterized ING Network Simulator

**Files:**
- Modify: `brian/chapter31.ipynb`
- Create: `tests/test_ch31_brian_ing_networks.py`

**Interfaces:**
- Consumes: `WB_EQS`, `ING_CONFIGS`, `build_ing_inputs(...)`, `solve_tau_dq(...)`.
- Produces: `simulate_ing_network(name, duration=500*b2.ms, dt=0.01*b2.ms, num_i=100, seed=None, record=False) -> dict`.
- Result keys: `t_i_ms`, `i_i`, `duration_ms`, `num_i`, `config`, and optionally `state`.

- [ ] **Step 1: Write failing structural and numeric tests**

```python
from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = ROOT / "python" / "31_ING_Rhythms"


@pytest.mark.parametrize("name", [f"ING_{index}" for index in range(1, 11)])
def test_brian_ing_network_matches_python_population_activity(name):
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    result = ns.simulate_ing_network(name)
    assert len(result["t_i_ms"]) > 500
    counts = np.bincount(result["i_i"].astype(int), minlength=result["num_i"])
    assert np.count_nonzero(counts) > 0.9 * result["num_i"]

    if name != "ING_1":
        py = load_python_port(PYTHON_BASE / name / "main.py")
        brian_rate = len(result["t_i_ms"]) / (result["num_i"] * result["duration_ms"] / 1000.0)
        python_rate = len(py.t_i_spikes) / (py.num_i * py.t_final / 1000.0)
        assert np.isclose(brian_rate, python_rate, rtol=0.20)


def test_ing_network_rejects_unknown_name():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    with pytest.raises(ValueError, match="ING_99"):
        ns.simulate_ing_network("ING_99", duration=20*ns.b2.ms, num_i=10)
```

- [ ] **Step 2: Run one representative test and verify RED**

Run: `python3 -m pytest tests/test_ch31_brian_ing_networks.py::test_brian_ing_network_matches_python_population_activity -k ING_2 -v`

Expected: FAIL because `simulate_ing_network` is absent.

- [ ] **Step 3: Implement chemical inhibition in Brian2**

Create a WB `NeuronGroup`, assign constants and generated initial states, then
map every nonzero `g_ii_ms[pre, post]` entry into a Brian synapse:

```python
pre, post = np.nonzero(inputs["g_ii_ms"])
chemical = b2.Synapses(
    cells,
    cells,
    model="g : siemens\nI_chem_post = g*s_pre*(-75*mV-v_post) : amp (summed)",
)
chemical.connect(i=pre, j=post)
chemical.g = inputs["g_ii_ms"][pre, post] * b2.msiemens
```

- [ ] **Step 4: Implement diffusive gap-junction coupling**

Map the symmetric matrix as directed Brian edges. The matrix is indexed as
`[post, pre]`, matching `G_gap @ v - row_sum*v` in the Python reference:

```python
post, pre = np.nonzero(inputs["g_gap_ms"])
gap = b2.Synapses(
    cells,
    cells,
    model="g : siemens\nI_gap_post = g*(v_pre-v_post) : amp (summed)",
)
gap.connect(i=pre, j=post)
gap.g = inputs["g_gap_ms"][post, pre] * b2.msiemens
```

When the relevant matrix is empty, still create the `Synapses` object but do
not connect edges. Use a `SpikeMonitor`; optionally add a `StateMonitor` for
`v`, `q`, and `s`; return plain NumPy values. Reject unknown configuration
names before calling `start_scope()`.

- [ ] **Step 5: Run the representative chemical and electrical cases**

Run: `python3 -m pytest tests/test_ch31_brian_ing_networks.py -k 'ING_2 or ING_6 or unknown' -v`

Expected: 3 passed.

- [ ] **Step 6: Run all ten network cases**

Run: `python3 -m pytest tests/test_ch31_brian_ing_networks.py -v`

Expected: 11 passed.

- [ ] **Step 7: Add guarded raster cells for all ten examples**

Define one plotting helper:

```python
def plot_ing_raster(result, ax=None, closeup=False):
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))
    ax.plot(result["t_i_ms"], result["i_i"], ".k", markersize=2)
    if closeup:
        ax.set_xlim(result["duration_ms"]-100.0, result["duration_ms"])
        ax.set_ylim(0, 21)
    else:
        ax.set_xlim(0, result["duration_ms"])
        ax.set_ylim(0, result["num_i"]+1)
    ax.set_xlabel("t [ms]")
    return ax
```

Under `if __name__ == "__main__"`, run `ING_1` through `ING_10`; pass
`closeup=True` for `ING_8`, `ING_9`, and `ING_10`.

- [ ] **Step 8: Commit Task 5 only**

```bash
git add brian/chapter31.ipynb tests/test_ch31_brian_ing_networks.py
git commit -m "feat: add Chapter 31 Brian ING network variants"
```

---

### Task 6: E/I Entrainment Networks

**Files:**
- Modify: `brian/chapter31.ipynb`
- Create: `tests/test_ch31_brian_ing_entrainment.py`

**Interfaces:**
- Consumes: `WB_EQS`, `solve_tau_dq(...)`, and the matrix orientation and normalization established by `build_ing_inputs(...)`.
- Produces: `RTM_EQS` and `EI_WB_EQS` equation strings.
- Produces: `simulate_ing_entrainment(variant, e_drive=None, duration=500*b2.ms, dt=0.01*b2.ms, num_e=400, num_i=100, seed=63806) -> dict`.
- Result keys: `t_e_ms`, `i_e`, `t_i_ms`, `i_i`, `duration_ms`, `num_e`, `num_i`, and `config`.

- [ ] **Step 1: Write failing entrainment tests**

```python
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_as_module, load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = ROOT / "python" / "31_ING_Rhythms"


def _active_fraction(indices, count):
    return np.count_nonzero(np.bincount(np.asarray(indices, dtype=int), minlength=count)) / count


def test_brian_ing_entraining_e_cells_matches_python_structure():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    py = load_python_port(PYTHON_BASE / "ING_ENTRAINING_E_CELLS" / "main.py")
    result = ns.simulate_ing_entrainment("ING_ENTRAINING_E_CELLS")
    assert _active_fraction(result["i_e"], result["num_e"]) > 0.9
    assert _active_fraction(result["i_i"], result["num_i"]) > 0.9
    assert np.isclose(len(result["t_e_ms"]), len(py.t_e_spikes), rtol=0.25)
    assert np.isclose(len(result["t_i_ms"]), len(py.t_i_spikes), rtol=0.25)


def test_brian_ing_entraining_e_cells_drive_sweep_has_three_results():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter31.ipynb")
    results = [
        ns.simulate_ing_entrainment("ING_ENTRAINING_E_CELLS_2", e_drive=drive)
        for drive in (1.9, 2.0, 2.1)
    ]
    assert len(results) == 3
    assert all(len(result["t_e_ms"]) > 500 for result in results)
    assert all(len(result["t_i_ms"]) > 200 for result in results)
```

- [ ] **Step 2: Run the test and verify RED**

Run: `python3 -m pytest tests/test_ch31_brian_ing_entrainment.py -v`

Expected: FAIL because `simulate_ing_entrainment` is absent.

- [ ] **Step 3: Add RTM equations and exact variant configurations**

Define the entrainment WB equation variant and the RTM equations:

```python
EI_WB_EQS = (
    WB_EQS
    .replace("i_ext + I_chem + I_gap", "i_ext + I_from_e + I_from_i + I_gap")
    .replace("I_chem : amp", "I_from_e : amp\nI_from_i : amp")
)

RTM_EQS = """
dv/dt = (g_l*(E_l-v) + g_k*n**4*(E_k-v)
         + g_na*m_inf**3*h*(E_na-v) + i_ext + I_from_e + I_from_i)/C : volt
dh/dt = (h_inf-h)/tau_h : 1
dn/dt = (n_inf-n)/tau_n : 1
dq/dt = 0.5*(1+tanh(v/(10*mV)))*(1-q)/(0.1*ms) - q/tau_dq : 1
ds/dt = q*(1-s)/tau_r - s/tau_d : 1
m_inf = alpha_m/(alpha_m+beta_m) : 1
h_inf = alpha_h/(alpha_h+beta_h) : 1
n_inf = alpha_n/(alpha_n+beta_n) : 1
alpha_m = 0.32/mV*(v+54*mV)/(1-exp(-(v+54*mV)/(4*mV)))/ms : Hz
beta_m = 0.28/mV*(v+27*mV)/(exp((v+27*mV)/(5*mV))-1)/ms : Hz
alpha_h = 0.128*exp(-(v+50*mV)/(18*mV))/ms : Hz
beta_h = 4/(1+exp(-(v+27*mV)/(5*mV)))/ms : Hz
alpha_n = 0.032/mV*(v+52*mV)/(1-exp(-(v+52*mV)/(5*mV)))/ms : Hz
beta_n = 0.5*exp(-(v+57*mV)/(40*mV))/ms : Hz
tau_h = 1/(alpha_h+beta_h) : second
tau_n = 1/(alpha_n+beta_n) : second
i_ext : amp
I_from_e : amp
I_from_i : amp
tau_dq : second
tau_r : second
tau_d : second
C : farad (constant)
g_l : siemens (constant)
g_k : siemens (constant)
g_na : siemens (constant)
E_l : volt (constant)
E_k : volt (constant)
E_na : volt (constant)
"""
```

Assign RTM constants `C=1 uF`, `g_l=0.1 mS`, `g_k=80 mS`,
`g_na=100 mS`, `E_l=-67 mV`, `E_k=-100 mV`, and `E_na=50 mV`.

The first variant uses:

```python
dict(num_e=400, num_i=100, sigma_e=0.10, sigma_i=0.05,
     mean_e=1.5, mean_i=1.5,
     g_hat_ee=0.0, g_hat_ei=0.0, g_hat_ie=0.50, g_hat_ii=0.50,
     p_ee=1.0, p_ei=0.5, p_ie=0.5, p_ii=0.5,
     g_hat_gap=0.1, p_gap=0.05)
```

The second uses zero heterogeneity, `mean_i=1.5`, selectable `mean_e` in
`(1.9, 2.0, 2.1)`, all four probabilities equal to `1.0`, chemical strengths
`(0.0, 0.0, 0.5, 0.5)`, and gap parameters `(0.1, 0.05)`.

- [ ] **Step 4: Implement the two-population simulator**

Create RTM and WB `NeuronGroup` objects. Generate and assign the four chemical
matrices with the Python normalization denominators:

- E-to-E: `g_hat_ee/(num_e*p_ee)`;
- E-to-I: `g_hat_ei/(num_e*p_ei)`;
- I-to-E: `g_hat_ie/(num_i*p_ie)`;
- I-to-I: `g_hat_ii/(num_i*p_ii)`.

Use continuous summed synaptic currents with `s_pre`, create the inhibitory
gap junctions with the following directed-edge mapping, initialize RTM and WB
populations from deterministic phase samples, and monitor spikes from both
groups:

```python
def add_continuous_synapses(source, target, weights_ms, reversal, target_variable):
    pre, post = np.nonzero(weights_ms)
    model = (
        "g : siemens\n"
        "E_rev : volt (constant)\n"
        f"{target_variable}_post = g*s_pre*(E_rev-v_post) : amp (summed)"
    )
    synapses = b2.Synapses(source, target, model=model)
    synapses.connect(i=pre, j=post)
    synapses.g = weights_ms[pre, post] * b2.msiemens
    synapses.E_rev = reversal
    return synapses


gap_post, gap_pre = np.nonzero(g_gap_ms)
gap = b2.Synapses(
    inhibitory,
    inhibitory,
    model="g : siemens\nI_gap_post = g*(v_pre-v_post) : amp (summed)",
)
gap.connect(i=gap_pre, j=gap_post)
gap.g = g_gap_ms[gap_post, gap_pre] * b2.msiemens
```

Use `add_continuous_synapses` four times: E-to-E targets `I_from_e` at
`0 mV`; E-to-I targets `I_from_e` at `0 mV`; I-to-E targets `I_from_i` at
`-75 mV`; and I-to-I targets `I_from_i` at `-75 mV`. Return only plain arrays
and scalar metadata.

- [ ] **Step 5: Run the entrainment tests**

Run: `python3 -m pytest tests/test_ch31_brian_ing_entrainment.py -v`

Expected: 2 passed.

- [ ] **Step 6: Add guarded raster plots**

For the first variant, draw I spikes in blue, E spikes shifted by `num_i` in
red, and the dashed population boundary. For the second, create three stacked
panels labeled with E drive `1.9`, `2.0`, and `2.1`, following the Python
figure limits and colors.

- [ ] **Step 7: Commit Task 6 only**

```bash
git add brian/chapter31.ipynb tests/test_ch31_brian_ing_entrainment.py
git commit -m "feat: add Chapter 31 Brian E-I entrainment"
```

---

### Task 7: Compatibility Status, Aggregate Verification, and Figures

**Files:**
- Modify: `PROGRESS.md`
- Modify: `brian/chapter31.ipynb`
- Test: `tests/test_ch31_brian_shared.py`
- Test: `tests/test_ch31_brian_1_cell_ing.py`
- Test: `tests/test_ch31_brian_1_cell_ing_condition_numbers.py`
- Test: `tests/test_ch31_brian_ing_configuration.py`
- Test: `tests/test_ch31_brian_ing_networks.py`
- Test: `tests/test_ch31_brian_ing_entrainment.py`

**Interfaces:**
- Consumes: every Chapter 31 Brian function and result contract.
- Produces: final Chapter 31 progress classification and visually checked notebook outputs.

- [ ] **Step 1: Add a notebook compatibility note**

Near the notebook title, add a markdown table listing all 16 Python examples.
Mark the 14 simulation examples `Brian2` and explain that the two abstract
pulse-coupling examples are analytical phase-map calculations available in
the Python and MATLAB directories, not Brian2 simulations.

- [ ] **Step 2: Run all Chapter 31 Brian tests**

Run:

```bash
python3 -m pytest -q \
  tests/test_ch31_brian_shared.py \
  tests/test_ch31_brian_1_cell_ing.py \
  tests/test_ch31_brian_1_cell_ing_condition_numbers.py \
  tests/test_ch31_brian_ing_configuration.py \
  tests/test_ch31_brian_ing_networks.py \
  tests/test_ch31_brian_ing_entrainment.py
```

Expected: all tests pass with zero failures.

- [ ] **Step 3: Run the existing Python-to-MATLAB Chapter 31 tests**

Run: `python3 -m pytest -q tests/test_ch31_*.py`

Expected: all Chapter 31 Python and Brian tests pass; MATLAB-only checks may
skip if MATLAB is unavailable.

- [ ] **Step 4: Execute the notebook at full duration**

Run:

```bash
jupyter nbconvert --to notebook --execute brian/chapter31.ipynb \
  --ExecutePreprocessor.timeout=1800 \
  --output chapter31.executed.ipynb \
  --output-dir /tmp
```

Expected: exit 0 with `/tmp/chapter31.executed.ipynb` created.

- [ ] **Step 5: Visually compare every generated figure**

Compare the single-cell trace, ten ING rastergrams, the first entrainment
raster, and the three-panel drive sweep against the corresponding PNGs under
`python/31_ING_Rhythms`. Check population ordering, time windows, active-cell
coverage, synchrony/desynchrony pattern, axis limits, and colors. If a visual
discrepancy reflects model behavior rather than styling, stop and fix the
model with a failing focused test before continuing.

- [ ] **Step 6: Update `PROGRESS.md`**

Set Brian status to `[x]` for:

```text
1_CELL_ING
1_CELL_ING_CONDITION_NUMBERS
ING_1
ING_2
ING_3
ING_4
ING_5
ING_6
ING_7
ING_8
ING_9
ING_10
ING_ENTRAINING_E_CELLS
ING_ENTRAINING_E_CELLS_2
```

Set Brian status to `n/a` for:

```text
ABSTRACT_PULSE_COUPLING_INH
ABSTRACT_PULSE_COUPLING_INH_2
```

- [ ] **Step 7: Run final repository checks**

Run:

```bash
python3 -m pytest -q tests/test_ch31_*.py
git diff --check
git status --short
```

Expected: tests pass, `git diff --check` is silent, and only Chapter 31 files
from this task plus pre-existing unrelated changes appear in status.

- [ ] **Step 8: Commit final Chapter 31 status**

```bash
git add brian/chapter31.ipynb PROGRESS.md
git commit -m "docs: complete Chapter 31 Brian coverage"
```
