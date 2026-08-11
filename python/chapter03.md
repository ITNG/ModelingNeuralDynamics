# The classical HH ODEs

## Overview

This chapter isolates the voltage-dependent gate functions used by the
classical HH model. Rather than simulating a trace, it plots the equilibrium
gate positions and their response speeds across voltage.

## Core ideas

For a gate, an opening rate and a closing rate determine both where it settles
and how quickly it gets there. Activation gates increase availability with
depolarization in their relevant range; inactivation can move in the opposite
direction. These curves explain why the same voltage equation can produce a
brief spike instead of unbounded depolarization.

## Essential model

For any gate \(x\), the implementation can be written as

\[
\frac{dx}{dt}=\frac{x_\infty(V)-x}{\tau_x(V)},\qquad
x_\infty(V)=\frac{\alpha_x(V)}{\alpha_x(V)+\beta_x(V)},\qquad
\tau_x(V)=\frac{1}{\alpha_x(V)+\beta_x(V)}.
\]

Here \(x\) is one of the dimensionless gates \(m\), \(h\), or \(n\); \(V\)
is membrane voltage; \(x_\infty\) is its steady-state value; \(\tau_x\) is
its voltage-dependent time constant; and \(\alpha_x\) and \(\beta_x\) are
the opening and closing rates, respectively.

## Code examples

The example now lives in [`chapter03.ipynb`](chapter03.ipynb):
`simulate_hh_gating_variables` evaluates `m_inf`, `h_inf`, `n_inf`, and
their time constants over a voltage grid, plotted as a six-panel figure.

## What to look for

Compare the left-column steady-state curves with the right-column time
constants. A high \(m_\infty\) does not mean immediate sodium activation:
the corresponding \(\tau_m\) says how quickly it can approach that level.

## Suggested order

1. Inspect the rate functions `alpha_*` and `beta_*`.
2. Relate each pair of rates to \(x_\infty\) and \(\tau_x\).
3. Return to the voltage equation in Chapter 01 and identify each gate's
   conductance exponent.

## Prerequisites and related chapters

Chapter 01 introduces the full HH balance equation. Chapter 04 uses the same
rates during numerical integration. Familiarity with exponential functions
and first-order relaxation is helpful.

## Running the examples

Open [`chapter03.ipynb`](chapter03.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom.
