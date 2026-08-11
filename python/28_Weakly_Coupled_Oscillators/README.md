# Weakly coupled oscillators

## Overview

When coupling is weak, the full neuron equations can be reduced to slowly
changing phases. These examples plot interaction and difference functions and
simulate identical and heterogeneous oscillator pairs to relate zeros of the
reduced equation to stable phase differences.

## Core ideas

Weak coupling preserves each oscillator's basic cycle while shifting its phase
slightly each period. The interaction function summarizes that shift. Subtracting
the two phase equations gives a difference function: its zeros are locked phase
differences, and its local sign or derivative determines which are stable.
Intrinsic frequency mismatch adds a drift term and can move or remove a locked
state.

## Essential model

For phases $\theta_1,\theta_2$, a weak-coupling reduction has the form

$$
\dot\theta_i=\omega_i+\varepsilon H(\theta_j-\theta_i).
$$

With $\psi=\theta_2-\theta_1$, the difference equation is
$\dot\psi=\omega_2-\omega_1+\varepsilon D(\psi)$. A zero with restoring
flow on both sides is a stable phase difference.

## Code examples

- [`PLOT_D_TWO_FIXED_POINTS`](PLOT_D_TWO_FIXED_POINTS/) plots a difference
  function with two phase-locked fixed points.
- [`WEAKLY_COUPLED_1`](WEAKLY_COUPLED_1/) simulates an identical weakly coupled
  pair for one interaction function.
- [`WEAKLY_COUPLED_2`](WEAKLY_COUPLED_2/) changes the interaction function to
  show a different phase-difference evolution.
- [`WEAKLY_COUPLED_HETEROGENEOUS_1`](WEAKLY_COUPLED_HETEROGENEOUS_1/) adds
  oscillator heterogeneity and shows its effect on locking.

## What to look for

In `PLOT_D_TWO_FIXED_POINTS`, identify each zero and inspect the direction of
the flow on either side before naming it stable. Compare the long-time phase
difference in the two identical-pair simulations with those predictions. In the
heterogeneous simulation, separate a genuine locked offset from ongoing phase
drift caused by frequency mismatch.

## Suggested order

1. Run `PLOT_D_TWO_FIXED_POINTS`.
2. Run `WEAKLY_COUPLED_1` and `WEAKLY_COUPLED_2`.
3. Run `WEAKLY_COUPLED_HETEROGENEOUS_1` and compare it with the identical
   cases.

## Prerequisites and related chapters

Chapter 25 derives PRCs and interaction functions, and Chapters 26--27 use
event maps with pulse coupling and delay. Chapter 29 gives a more detailed
view of synchronous-state stability under inhibitory perturbations.

## Running the examples

Run `python main.py` from the desired immediate example directory. These are
small NumPy and Matplotlib simulations that save their figures locally.
