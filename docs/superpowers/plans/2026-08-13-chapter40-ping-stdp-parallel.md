# Chapter 40 PING-with-STDP Parallel Acceleration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Accelerate `simulate_ping_with_stdp()` through exact row-parallel STDP updates and compiled population initialization.

**Architecture:** Preserve the sequential Heun timestep loop while parallelizing independent synapse rows inside each STDP derivative evaluation. Compile initialization and scope the Numba worker cap to the full-network loop call.

**Tech Stack:** Python, NumPy, Numba, pytest, Jupyter notebook JSON

## Global Constraints

- Preserve `simulate_ping_with_stdp` parameters and return structure.
- Preserve `dt=0.01`, both Heun stages, random-number order, and all STDP equations.
- Restore the previous Numba thread setting after simulation, including when the loop raises.

---

### Task 1: Compile initialization and parallelize the STDP kernel

**Files:**
- Modify: `tests/test_ch40_stdp.py`
- Modify: `python/chapter40.ipynb`

**Interfaces:**
- Consumes: `rtm_init_population(i_ext, phi_vec)`, `_g_ee_derivative_s(...)`, and `simulate_ping_with_stdp(...)`
- Produces: unchanged numerical results through Numba nopython and parallel specializations

- [x] **Step 1: Write failing compilation and equivalence tests**

Add a short real simulation test that requires nopython initialization and a parallel STDP kernel, plus a direct small-matrix kernel test with hand-computed serial equations.

- [x] **Step 2: Verify the tests fail for missing compilation modes**

Run the two targeted tests and confirm failure because initialization is Python and the STDP kernel has no parallel overload.

- [x] **Step 3: Implement the minimal exact optimization**

Import `prange`, `get_num_threads`, and `set_num_threads`; compile `rtm_init_population`; use explicit Numba dtypes; change only the dense row loop in `_g_ee_derivative_s` to `prange`; enable `parallel=True`; and restore the thread setting in `finally`.

- [x] **Step 4: Run targeted and Chapter 40 tests**

Run both targeted tests followed by `MPLCONFIGDIR=/tmp python3 -m pytest tests/test_ch40_stdp.py -v`.

- [x] **Step 5: Verify deterministic equivalence and benchmark**

Compare serial/parallel kernel values, compare representative simulation results, benchmark compile and warm default runs, validate notebook JSON, and run `git diff --check`.
