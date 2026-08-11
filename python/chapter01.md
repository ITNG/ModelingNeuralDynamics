# Modeling a single neuron

## Overview

This directory is the repository's first executable modeling example. It is
not a mirror of the book's vocabulary-only Chapter 1: it starts immediately
with a conductance-based Hodgkin--Huxley (HH) simulation and a voltage trace.

## Core ideas

A membrane stores charge, so its voltage changes only when the injected and
ionic currents do not balance. Sodium activation provides rapid positive
feedback; sodium inactivation and potassium activation provide delayed
negative feedback. The three gating variables make those conductances depend
on both voltage and recent history.

## Essential model

The voltage balance implemented here is

$$
C\frac{dV}{dt}=I_{\mathrm{ext}}-g_{\mathrm{Na}}m^3h(V-E_{\mathrm{Na}})
-g_{\mathrm{K}}n^4(V-E_{\mathrm{K}})-g_{\mathrm{L}}(V-E_{\mathrm{L}}).
$$

Here $V$ is membrane voltage, $t$ is time, $C$ is membrane
capacitance, $I_{\mathrm{ext}}$ is applied current, $g_{\mathrm{Na}}$,
$g_{\mathrm{K}}$, and $g_{\mathrm{L}}$ are maximal sodium, potassium,
and leak conductances, and $E_{\mathrm{Na}}$, $E_{\mathrm{K}}$, and
$E_{\mathrm{L}}$ are their reversal potentials. $m$, $h$, and $n$
are dimensionless sodium-activation, sodium-inactivation, and
potassium-activation gates. Each gate $x\in\{m,h,n\}$ follows
$dx/dt=\alpha_x(V)(1-x)-\beta_x(V)x$, where $\alpha_x$ and $\beta_x$
are voltage-dependent opening and closing rates.

## Code examples

The example now lives in [`chapter01.ipynb`](chapter01.ipynb):
`simulate_hh_voltage_trace` integrates the four HH state variables from
their voltage-dependent initial values, and an `ipywidgets` slider lets you
adjust `i_ext` interactively.

## What to look for

The trace is the compact consequence of current balance: a depolarizing input
first recruits $m$, then $h$ falls and $n$ rises to end the spike. Try
changing `i_ext` only after confirming how the same initial state is built
from `m_inf`, `h_inf`, and `n_inf`.

## Suggested order

1. Run `HH_VOLTAGE_TRACE` and identify the resting level, upstroke, and
   recovery.
2. Read the `derivative` function alongside the equation above.
3. Continue to the gating curves in Chapter 03 before changing rate laws.

## Prerequisites and related chapters

Comfort with first-order ODEs, units of current and voltage, and Python arrays
is enough. Chapter 03 separates the gates into steady-state and time-constant
curves; Chapter 04 focuses on numerical trajectories.

## Running the examples

Open [`chapter01.ipynb`](chapter01.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook. Run all cells top to bottom; the static
figure reproduces the book's plot, and the `interact(...)` cell below it
lets you adjust `i_ext` with a slider.
