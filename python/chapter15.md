# Canard explosions

## Overview

These examples resolve the abrupt growth of an oscillation in a slow-fast
system. FitzHugh--Nagumo parameter scans show the macro and micro views of a
canard explosion, while reduced HH and adaptation examples show related
multiple time-scale behavior.

## Core ideas

Near a canard, a small parameter change can move a trajectory from a tiny
oscillation to a large relaxation cycle. The trajectory follows an otherwise
repelling slow-manifold branch for an anomalously long interval. Macro scans
show the global branches; micro scans reveal the extremely narrow parameter
window where cycle amplitude changes sharply.

## Essential model

The FitzHugh--Nagumo system is

$$
\dot v=v-v^3/3-n+I,\qquad \dot n=(av-n)/\tau_n.
$$

The slow recovery variable $n$ has time scale $\tau_n$, $v$ is fast,
and $I$ is the scanned control parameter. In the adaptation example, an
additional slowly decaying feedback current shifts the effective drive during
the oscillation.

## Code examples

All six examples now live in one notebook, [`chapter15.ipynb`](chapter15.ipynb):
`simulate_canard` finds currents for selected FitzHugh-Nagumo amplitudes
and plots phase-plane and time-trace views; `simulate_canard_2` fixes a
current inside the narrow transition and plots its phase-plane orbit;
`simulate_fitzhugh_nagumo_macro` scans a broad current range for
equilibrium and cycle envelopes; `simulate_fitzhugh_nagumo_micro` resolves
the narrow near-critical scan and cycle amplitudes; `simulate_hh_reduced_bif_diag`
plots reduced HH fixed-point and stable/unstable cycle envelopes;
`simulate_mmos` adds slow adaptation to FitzHugh-Nagumo and plots its
mixed-mode voltage trace (with an `ipywidgets` slider). Several of these
cells are computationally expensive (minutes, not seconds) -- the notebook
notes this inline.

## What to look for

Compare the scale of `FITZHUGH_NAGUMO_MACRO` and
`FITZHUGH_NAGUMO_MICRO`: the latter is needed to see the amplitude jump.
`CANARD` makes the same transition visible as selected trajectories. In
`HH_REDUCED_BIF_DIAG`, distinguish the stable cycle envelope from the
backward-traced unstable one.

## Suggested order

1. Run `FITZHUGH_NAGUMO_MACRO`, then `FITZHUGH_NAGUMO_MICRO`.
2. Run `CANARD` and `CANARD_2` to inspect individual trajectories.
3. Continue with `HH_REDUCED_BIF_DIAG` and `MMOS`.

## Prerequisites and related chapters

Chapter 10 introduces slow-fast phase planes and Chapter 13 introduces Hopf
and unstable cycles. Chapter 14 supplies the reduced HH type-2 setting;
Chapter 19 uses a still slower current to make bursts.

## Running the examples

Open [`chapter15.ipynb`](chapter15.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom. The
canard-scan cells are noticeably slower than the single-trajectory ones.
