# Brian Chapter Coverage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Execute every Python notebook and supply validated representative Brian2
examples for all 23 compatible missing chapters.

**Architecture:** Shared mnd.brian neuron/dynamics/phase/network/plasticity modules;
short notebooks describe chosen experiments and plot results. Independent Python
chapter definitions serve as numerical references in tests only.

**Tech Stack:** existing uv environment, Brian2, NumPy/SciPy, nbclient, pytest.

**Spec:** docs/superpowers/specs/2026-10-28-brian-chapter-coverage-design.md

## Global Constraints

- Representative coverage, not full parity; 03 excluded as entirely static.
- Preserve scientific equations and state initialization; disclose deviations.
- No production AST loading or execution of Python reference notebooks.
- Empty tracked notebook metadata; source checkout installation for new helpers.
- Three cycles; no later questions or unsolicited commit/push.

## Review Focus

- Singular HH/RTM rates at exact voltages must remain finite (use exprel).
- Independent clocks/namespaces and repeated simulations must not share state.
- Spike timing, phase reset order, and transmission delays must match conventions.
- Plasticity gates/weights must be finite, bounded, and respond to activity.
- Missing/unequal population sizes must not reverse synaptic normalization.

### Task 1 — Baseline Python notebooks (python_validation)
- [x] Execute all 38 in fresh kernels with original defaults; retain outputs/logs.
- [x] Reproduce and fix any execution errors at source, run relevant regressions.
- [x] Rerun changed notebooks and report exact coverage and remaining limitations.

### Task 2 — Shared neurons and phase chapters (root)
Files: mnd/brian/models.py, mnd/brian/phases.py; tests/test_brian_models.py,
tests/test_brian_phases.py; brian/chapter25.ipynb through chapter29.ipynb.
- [x] Test factory fields, singular rates, integration vs independent SciPy, isolation.
- [x] Implement make_neurons(model,n=1,current=0,*,dt=.01,v0=None,
  instantaneous_m=None,parameters=None,extra_equations='',synaptic_current=None,
  namespace=None), numeric v/current, ODE /ms, explicit local clock/namespace.
- [x] Test phase oscillators against exact event references, weak phase ODE against
  Python trajectories, and LIF timing against closed form.
- [x] Implement shared perturbation/delayed-phase/weak-coupling/inhibition protocols;
  create five documented notebooks and execute them fully.

### Task 3 — Early dynamics chapters (brian_scope)
Files: mnd/brian/dynamics.py; tests/test_brian_dynamics.py; chapters10–18,22.
- [x] Numerical tests first against SciPy/Python reference dynamics and known roots.
- [x] Implement reusable FHN, normal forms, reduced neurons, Wilson–Cowan and theta.
- [x] Create ten representative notebooks, including FI and bistability experiments.
- [x] Execute every new notebook; retain figures and numerical validation artifacts.

### Task 4 — Network/drive/plasticity chapters (brian_networks)
Files: mnd/brian/networks.py, drives.py, plasticity.py as needed;
tests/test_brian_network_extensions.py, test_brian_plasticity.py;
chapters32–33,35–40.
- [x] Numerical tests first for source equations, projections, gates and drive means.
- [x] Implement shared smooth synapses, drives, adaptation and plasticity protocols.
- [x] Create eight representative notebooks with explicit coverage/omission notes.
- [x] Execute all notebooks and compare selected source trajectories/rates.

### Task 5 — Integration, documentation and final cycle (root)
- [x] Document complete Brian coverage/exclusions, setup, shared modules and examples.
- [x] Run all Brian notebooks with final modules; all Python notebooks verified.
- [x] Run default and explicit new/existing Brian regression suites; inspect failures.
- [x] Independent review, resolve material issues with tests, record final evidence.

## Progress ledger

- Baseline: HEAD3acfaab, only untracked .clinerules; existing19Brian notebooks and
  65Brian/241default tests previously verified. Python38 execution started.
- Scope:23 new notebooks; only03 excluded. User confirmed representativeexamples.
- Ownership:root models/phases/docs; earlyagent dynamics+10notebooks; networkagent
  networks/drives/plasticity+8notebooks; Pythonagent only Python source/artifacts.

- Task 2 complete: 12 neuron/phase tests passed; exact Erisir singularity
  voltages subsequently verified in 7 model tests; all 5 phase notebooks passed.
- Task 3 complete: 8 numerical tests and 10 fresh notebooks passed.
- Task 4 complete: 9 initial tests and 8 fresh notebooks passed. Review added
  actual Brian Synapses derivative and facilitation trajectory checks; all
  6 plasticity tests passed, including the two independent checks.
- Final Brian execution: 42/42 passed; 181 nonempty code cells, 102 figures;
  no errors or numerical overflow/invalid warnings; source hashes match.
- Python execution: 38/38 passed against final source; 616/616 code cells,
  321 figures, zero error outputs and zero stderr. Changed chapters 13, 17 and
  34 were fully rerun with original defaults. Chapter 34's reversed E-cell
  synaptic conductance signs were fixed against the MATLAB reference.
- Final default suite: 248 passed, 3 MATLAB-dependent skips, after the chapter
  34 correction and its regression. Ruff and locked uv sync pass.
- Independent review found no material defects in the new Brian equations;
  both suggested validation gaps now have direct numerical checks.

- Final Brian suite: 94 passed; the two additional review tests passed in the
  6-test plasticity run. Dynamics tests now restore their backend setting;
  verified by an 8-test in-process run and a before/after preference assertion.
- Complete: all 38 Python and 42 Brian notebooks verified. No commit or push.
