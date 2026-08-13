# Weakly coupled oscillators

## Overview

When coupling is weak, the full neuron equations can be reduced to slowly
changing phases. The examples plot interaction and difference functions and
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

All four examples now live in one notebook, [`chapter28.ipynb`](chapter28.ipynb):
`find_two_fixed_points` builds the difference function `D_two_fixed_points`
(from the single-pulse-like PRC `g_0`) and locates its unstable and stable
zeros above $\psi=0.5$ by bisection; `plot_d_two_fixed_points` draws it.
`simulate_weakly_coupled_event_driven` and `simulate_weakly_coupled_de` are
the shared machinery for an identical weakly coupled pair -- the first
advances each oscillator exactly to its next spike (`ceiling` handles the
"just spiked" edge case) and applies the weak-coupling kick from a given `g`
at each spike, the second integrates the reduced difference equation with
Heun's method. `simulate_weakly_coupled_1` (interaction function `wc1_g`,
$\varphi^2(1-\varphi)$) and `simulate_weakly_coupled_2` (interaction function
`wc2_g`, $\varphi(1-\varphi)^3$) each call both for a chosen `epsilon`, and
`plot_weakly_coupled` overlays the event-driven and reduced-equation traces.
`simulate_weakly_coupled_heterogeneous_1` adds a period mismatch (B's period
is $T_B=1+\varepsilon c$) via `simulate_heterogeneous_event_driven` and
`simulate_heterogeneous_de`, and `plot_weakly_coupled_heterogeneous` compares
two detunings `c`.

## What to look for

In `find_two_fixed_points`, identify each zero and inspect the direction of
the flow on either side before naming it stable. Compare the long-time phase
difference in the two identical-pair simulations (`simulate_weakly_coupled_1`,
`simulate_weakly_coupled_2`) with those predictions. In the heterogeneous
simulation, separate a genuine locked offset from ongoing phase drift caused
by frequency mismatch.

## Suggested order

1. Run `find_two_fixed_points` / `plot_d_two_fixed_points`.
2. Run `simulate_weakly_coupled_1` and `simulate_weakly_coupled_2`.
3. Run `simulate_weakly_coupled_heterogeneous_1` and compare it with the
   identical cases.

## Prerequisites and related chapters

Chapter 25 derives PRCs and interaction functions, and Chapters 26-27 use
event maps with pulse coupling and delay. Chapter 29 gives a more detailed
view of synchronous-state stability under inhibitory perturbations.

## Running the examples

Open [`chapter28.ipynb`](chapter28.ipynb) in Jupyter, or via the Colab badge
at the top of the notebook, and run all cells top to bottom. These are small
NumPy/Matplotlib simulations that display their figures inline.
