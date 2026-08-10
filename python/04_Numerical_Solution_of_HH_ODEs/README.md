# Numerical solution of HH ODEs

## Overview

These examples use numerical integration to turn the HH equations into
time-series and state-space plots. Together they show a transient approaching
a repeating orbit and how the evolving gates produce refractoriness.

## Core ideas

An ODE solver approximates a continuous trajectory at a finite set of times.
After a transient, a periodically spiking HH state can approach a limit cycle.
During and after a spike, sodium inactivation and potassium activation leave
the state temporarily less excitable; a pulse therefore has a different
effect at different phases.

## Essential model

The solver advances the state vector

\[
\mathbf{x}=(V,m,n,h),\qquad \frac{d\mathbf{x}}{dt}=\mathbf{F}(\mathbf{x};I_{\mathrm{ext}}).
\]

Here \(V\) is membrane voltage; \(m\), \(n\), and \(h\) are the
dimensionless HH gates; \(t\) is time; \(I_{\mathrm{ext}}\) is applied
current; and \(\mathbf{F}\) is the HH voltage-and-gate right-hand side
defined in each `derivative` function. A limit cycle is a closed trajectory
in this state space that repeats after one period.

## Code examples

- [`HH_SOLUTION`](HH_SOLUTION/) integrates voltage and all three gates at a
  constant input and saves the two-panel result as `fig_4_1.png`.
- [`HH_LIMIT_CYCLE`](HH_LIMIT_CYCLE/) starts away from equilibrium, plots
  \(n\) against \(V\), and saves the projected cycle as `fig_4_2.png`.
- [`HH_REFRACTORINESS`](HH_REFRACTORINESS/) adds a short current pulse at
  several onset times and saves the comparison as `fig_4_1.png`.

## What to look for

`HH_SOLUTION` shows the initial transient before the repeating regime.
`HH_LIMIT_CYCLE` removes time from view so the returning state is visible.
In `HH_REFRACTORINESS`, compare the response to each red pulse marker: the
same perturbation can fail, delay, or advance a spike depending on gate state.

## Suggested order

1. Run `HH_SOLUTION` and identify its transient and later periodic segment.
2. Run `HH_LIMIT_CYCLE` to see that periodic segment as a closed projection.
3. Run `HH_REFRACTORINESS` to connect trajectory phase with excitability.

## Prerequisites and related chapters

Chapter 01 supplies the HH current balance and Chapter 03 explains its gate
functions. The LIF comparison in Chapter 07 offers a reset-based alternative
with no explicit ionic gates.

## Running the examples

Each directory has a `main.py` entry point. For example:

```bash
cd HH_SOLUTION
python main.py
```

The examples require SciPy, NumPy, and Matplotlib/Pylab and save figures in
their own working directories.
