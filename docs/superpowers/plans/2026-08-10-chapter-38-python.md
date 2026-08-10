# Chapter 38 Python Port Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port all five MATLAB Chapter 38 examples to readable, self-contained Python programs with reproducible figures and automated behavioral verification.

**Architecture:** Create `python/38_Gamma_Coherence/` with one folder per MATLAB example. Keep each example runnable through its own `main.py`, while placing the Hodgkin-Huxley gating equations and midpoint stepping used by three or more examples in a chapter-local `model.py`; this avoids premature changes to the public `mnd` package while keeping the five ports consistent. Deterministic coherence examples expose voltage, synaptic-gate, and spike arrays; stochastic PING examples expose seeded spike summaries and phase-response counts.

**Tech Stack:** Python 3, NumPy, Matplotlib, pytest, MATLAB R2020a reference scripts

## Global Constraints

- Mirror each MATLAB leaf folder as a Python leaf folder.
- Keep the code minimal and readable for beginners.
- Use the MATLAB midpoint method and `dt = 0.01` ms because stochastic event timing and spike thresholds are sensitive to the integration scheme.
- Seed NumPy with `63806`; because MATLAB and NumPy random streams differ, compare stochastic examples by invariant/statistical outputs rather than exact spike times.
- A sub-example is complete only after its figure is regenerated and an automated numeric or behavioral comparison passes.
- Do not modify or stage unrelated existing changes on `develop_python`.

---

### Task 1: Shared Chapter 38 Neuron and Synapse Model

**Files:**
- Create: `python/38_Gamma_Coherence/model.py`
- Test: `tests/test_ch38_gamma_coherence.py`

**Interfaces:**
- Produces: `m_e_inf(v)`, `h_e_inf(v)`, `tau_h_e(v)`, `n_e_inf(v)`, `tau_n_e(v)`, the corresponding `_i` functions, `tau_d_q(tau_d, tau_r, tau_peak)`, `rtm_initial_state(i_ext, phases)`, `wb_initial_state(i_ext, phases)`, and `detect_downward_spikes(v_old, v_new, t_old, t_new)`.

- [ ] **Step 1: Write the failing shared-model test**

```python
def test_shared_model_matches_matlab_gate_values(ch38_model):
    v = np.array([-75.0, -60.0, -20.0])
    np.testing.assert_allclose(ch38_model.m_e_inf(v), [0.00263048, 0.05624786, 0.94432979], rtol=1e-5)
    np.testing.assert_allclose(ch38_model.m_i_inf(v), [0.00804324, 0.05293249, 0.81665924], rtol=1e-5)
    assert ch38_model.tau_d_q(3.0, 0.3, 0.3) > 0
```

- [ ] **Step 2: Run the test and verify it fails because `model.py` is absent**

Run: `python3 -m pytest tests/test_ch38_gamma_coherence.py::test_shared_model_matches_matlab_gate_values -v`

- [ ] **Step 3: Implement the shared equations and initial-state helpers**

Translate the equations in `matlab/38/*/{m,h,n,tau}_*.m`, `tau_d_q_function.m`, `rtm_init.m`, and `wb_init.m`. Functions must accept scalars or NumPy arrays and return NumPy-compatible values.

- [ ] **Step 4: Run the focused test and verify it passes**

Run: `python3 -m pytest tests/test_ch38_gamma_coherence.py::test_shared_model_matches_matlab_gate_values -v`

### Task 2: Deterministic Gamma Coherence Examples

**Files:**
- Create: `python/38_Gamma_Coherence/GAMMA_COHERENCE_1/main.py`
- Create: `python/38_Gamma_Coherence/GAMMA_COHERENCE_2/main.py`
- Modify: `tests/test_ch38_gamma_coherence.py`

**Interfaces:**
- Produces: `simulate(...) -> dict[str, np.ndarray | float]` in each module.
- GAMMA_COHERENCE_1 result keys: `t`, `v_e`, `v_i`, `s_i`, `t_e_spikes`, `t_i_spikes`, `i_main`, `i_dist`.
- GAMMA_COHERENCE_2 result keys: `coupled`, `mean_inhibition`, `mean_s_i`, where each run contains the Task 1 trace keys.

- [ ] **Step 1: Add failing tests for deterministic trace shape and MATLAB spike summaries**

```python
def test_gamma_coherence_1_matches_matlab_summary():
    result = load_ch38("GAMMA_COHERENCE_1").simulate(t_final=200.0)
    assert result["t"].shape == result["v_e"].shape
    np.testing.assert_allclose(result["t_e_spikes"], MATLAB_GAMMA_1_E_SPIKES, rtol=1e-2, atol=0.05)

def test_gamma_coherence_2_uses_mean_inhibition():
    result = load_ch38("GAMMA_COHERENCE_2").simulate(t_final=200.0)
    assert 0.0 < result["mean_s_i"] < 1.0
    assert len(result["coupled"]["t_e_spikes"]) != len(result["mean_inhibition"]["t_e_spikes"])
```

- [ ] **Step 2: Run both tests and verify they fail because the ports are absent**

Run: `python3 -m pytest tests/test_ch38_gamma_coherence.py -k 'gamma_coherence_1 or gamma_coherence_2' -v`

- [ ] **Step 3: Implement the two scalar E-I midpoint simulations and plotting entry points**

Use the exact MATLAB parameters and periodic-input normalization. Store traces needed for numeric checks and recreate each source `figure.pdf` as `fig.png` when run as a script.

- [ ] **Step 4: Generate MATLAB reference summaries and copy the numeric arrays into the tests**

Run each original with a temporary MATLAB command that saves spike times and `mean_s_i` to CSV under `/tmp`; do not modify tracked MATLAB sources.

- [ ] **Step 5: Verify both deterministic examples and generate their figures**

Run: `MPLBACKEND=Agg python3 -m pytest tests/test_ch38_gamma_coherence.py -k 'gamma_coherence_1 or gamma_coherence_2' -v`

Run: `cd python/38_Gamma_Coherence/GAMMA_COHERENCE_1 && MPLBACKEND=Agg python3 main.py`

Run: `cd python/38_Gamma_Coherence/GAMMA_COHERENCE_2 && MPLBACKEND=Agg python3 main.py`

### Task 3: Seeded Poisson PING Simulation

**Files:**
- Modify: `python/38_Gamma_Coherence/model.py`
- Create: `python/38_Gamma_Coherence/POISSON_PING_3_PLUS_GREEN/main.py`
- Modify: `tests/test_ch38_gamma_coherence.py`

**Interfaces:**
- Produces: `simulate(seed=63806, t_final=500.0) -> dict` with `t_e_spikes`, `i_e_spikes`, `t_i_spikes`, `i_i_spikes`, `f_hat_e`, and `f_hat_i`.
- Produces: a reusable chapter-local population midpoint kernel accepting an optional deterministic pulse callable and seeded RNG.

- [ ] **Step 1: Add a failing reproducibility and activity test**

```python
def test_poisson_ping_is_reproducible_and_active():
    module = load_ch38("POISSON_PING_3_PLUS_GREEN")
    first = module.simulate(seed=63806, t_final=100.0)
    second = module.simulate(seed=63806, t_final=100.0)
    np.testing.assert_array_equal(first["t_e_spikes"], second["t_e_spikes"])
    assert first["f_hat_e"] > 0
    assert first["f_hat_i"] > 0
```

- [ ] **Step 2: Run it and verify it fails because the simulation is absent**

Run: `python3 -m pytest tests/test_ch38_gamma_coherence.py::test_poisson_ping_is_reproducible_and_active -v`

- [ ] **Step 3: Implement the all-to-all Poisson PING kernel and green raster example**

Exploit the MATLAB examples' `p = 1` and uniform weights by using population sums instead of allocating four dense connectivity matrices. Preserve per-neuron Poisson arrivals and seeded initial phases.

- [ ] **Step 4: Compare full-run firing rates and raster structure with MATLAB**

Run MATLAB headlessly to record `f_hat_e`, `f_hat_i`, and spike counts. Assert broad statistical tolerances that accommodate different random generators while still detecting broken dynamics.

- [ ] **Step 5: Run the test and generate the figure**

Run: `MPLBACKEND=Agg python3 -m pytest tests/test_ch38_gamma_coherence.py::test_poisson_ping_is_reproducible_and_active -v`

Run: `cd python/38_Gamma_Coherence/POISSON_PING_3_PLUS_GREEN && MPLBACKEND=Agg python3 main.py`

### Task 4: Phase-Response Poisson PING Variants

**Files:**
- Create: `python/38_Gamma_Coherence/POISSON_PING_3_PLUS_PULSES/main.py`
- Create: `python/38_Gamma_Coherence/POISSON_PING_3_MISMATCHED_PULSES/main.py`
- Modify: `tests/test_ch38_gamma_coherence.py`

**Interfaces:**
- Each module produces `simulate_phases(seed=63806, phases=np.arange(0.1, 1.0, 0.1), t_final=500.0) -> dict` with `phases`, `spike_counts`, `fit_slope`, and `fit_intercept`.
- `PLUS_PULSES` uses `P = 31` ms; `MISMATCHED_PULSES` uses `P = 29` ms.

- [ ] **Step 1: Add failing phase-response tests**

```python
@pytest.mark.parametrize((name, period), [
    ("POISSON_PING_3_PLUS_PULSES", 31.0),
    ("POISSON_PING_3_MISMATCHED_PULSES", 29.0),
])
def test_poisson_ping_phase_response(name, period):
    module = load_ch38(name)
    result = module.simulate_phases(seed=63806, phases=np.array([0.1, 0.5, 0.9]), t_final=100.0)
    assert module.P == period
    assert result["spike_counts"].shape == (3,)
    assert np.all(result["spike_counts"] >= 0)
```

- [ ] **Step 2: Run the parameterized test and verify it fails because both modules are absent**

Run: `python3 -m pytest tests/test_ch38_gamma_coherence.py::test_poisson_ping_phase_response -v`

- [ ] **Step 3: Implement both thin variants over the shared population kernel**

Each phase must start from a deterministic child seed so results are reproducible and independent of the order in which phases are evaluated. Count spikes only from E-cells 1-5, matching MATLAB.

- [ ] **Step 4: Compare trends with MATLAB and render both phase-response figures**

Record MATLAB's nine spike counts and fitted slope for each period. Require the Python result to reproduce the qualitative slope sign and remain in the MATLAB count range with a deliberately broad tolerance.

- [ ] **Step 5: Run focused and full Chapter 38 tests**

Run: `MPLBACKEND=Agg python3 -m pytest tests/test_ch38_gamma_coherence.py -v`

### Task 5: Progress Tracking and Final Verification

**Files:**
- Modify: `PROGRESS.md`

**Interfaces:**
- Produces: Chapter 38 Python status of 5/5 complete; Brian status remains unchanged.

- [ ] **Step 1: Visually inspect all five generated figures against the MATLAB PDFs**

Check raster colors, axes, time spans, phase-response direction, and the coherence comparison panels.

- [ ] **Step 2: Update only the five Chapter 38 Python cells and the Python total**

Change each Chapter 38 Python `[ ]` to `[x]` and increment the Python total from `191/256` to `196/256`; leave Brian at `41/256`.

- [ ] **Step 3: Run the complete test suite**

Run: `MPLBACKEND=Agg python3 -m pytest -q`

- [ ] **Step 4: Audit the diff for unrelated files**

Run: `git status --short` and `git diff -- PROGRESS.md tests/test_ch38_gamma_coherence.py python/38_Gamma_Coherence docs/superpowers/plans/2026-08-10-chapter-38-python.md`
