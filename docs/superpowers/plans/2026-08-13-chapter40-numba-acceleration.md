# Chapter 40 Numba Acceleration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Accelerate the Chapter 40 three-cell PING simulations and E-to-E sweep by compiling the ODE right-hand side with Numba.

**Architecture:** Preserve `odeint` and all public notebook functions. Compile only the hot numerical channel-rate helpers and network derivative, using Numba-compatible matrix-vector operators.

**Tech Stack:** Python, NumPy, SciPy `odeint`, Numba, pytest, Jupyter notebook JSON

## Global Constraints

- Preserve `simulate_three_cell_ping` and `sweep_three_cell_ping_ee` signatures and return structures.
- Preserve adaptive `odeint` integration.
- Do not enable `fastmath`; numerical changes must remain within integration tolerances.

---

### Task 1: Protect and compile the three-cell PING hot path

**Files:**
- Modify: `tests/test_ch40_stdp.py`
- Modify: `python/chapter40.ipynb`

**Interfaces:**
- Consumes: `simulate_three_cell_ping(g_ee, ..., t_final, dt)` and `derivative_three_cell_ping`
- Produces: the same interfaces, with `derivative_three_cell_ping` compiled in Numba nopython mode

- [x] **Step 1: Write the failing test**

Add a test that runs a short simulation, checks its output shape and finite values, and asserts that `derivative_three_cell_ping.nopython_signatures` is populated.

- [x] **Step 2: Run the test to verify it fails**

Run: `MPLCONFIGDIR=/tmp python3 -m pytest tests/test_ch40_stdp.py::test_three_cell_ping_uses_compiled_rhs -v`

Expected: fail because the current Python function has no `nopython_signatures` attribute.

- [x] **Step 3: Write the minimal implementation**

Decorate the ten E/I channel-rate functions and `derivative_three_cell_ping` with `@njit`. Replace its four `np.matmul(matrix, vector)` calls with `matrix @ vector`. Update notebook prose that currently says this path does not use Numba.

- [x] **Step 4: Run targeted and Chapter 40 tests**

Run the targeted test, then `MPLCONFIGDIR=/tmp python3 -m pytest tests/test_ch40_stdp.py -v`.

- [x] **Step 5: Verify numerical equivalence and performance**

Compare representative weak/strong compiled outputs against an undecorated copy of the derivative and report maximum state and spike-time differences. Benchmark warm weak/strong calls and a short sweep.
