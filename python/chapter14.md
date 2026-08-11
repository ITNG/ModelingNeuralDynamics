# Model neurons of bifurcation type 2

## Overview

These examples study type-2 excitability, where firing begins at a nonzero
frequency through Hopf-related geometry. Reduced HH and Erisir models expose
fixed points, eigenvalues, attracting cycles, and repelling cycle boundaries.

## Core ideas

Unlike type-1 onset, type-2 onset has a finite firing frequency at threshold.
The reduced HH construction takes sodium activation to equilibrium and ties
inactivation to potassium activation. Eigenvalues show when an equilibrium
changes stability, while a repelling cycle can divide resting and spiking
basins.

## Essential model

The reduced HH equations are

$$
C\dot v=I_{\rm Na}(v,n)+I_{\rm K}(v,n)+I_L(v)+I,\qquad
\dot n=\alpha_n(v)(1-n)-\beta_n(v)n,
$$

with $m=m_\infty(v)$ and $h=0.83-n$. The Erisir reduction has the same
form but uses its rodent-neuron conductances and $h=0.36-n$. Here $v$ is
voltage, $n$ is a recovery gate, and $I$ is applied current.

## Code examples

All seven examples now live in one notebook, [`chapter14.ipynb`](chapter14.ipynb):
`simulate_erisir_reduced` compares three- and two-dimensional Erisir voltage
traces; `simulate_erisir_2d_fp` classifies Erisir reduced fixed points
across current; `simulate_hh_reduced_count_fp` reports the minimum and
maximum number of reduced HH fixed points over its current scan;
`simulate_hh_reduced_fixed_points` plots stable and unstable reduced HH
fixed-point branches; `simulate_hh_reduced_fp_evs` plots the real and
imaginary parts of fixed-point eigenvalues; `simulate_hh_reduced_repelling_cycle`
traces attracting and backward-integrated repelling cycles;
`simulate_hh_reduced_cycle_distance` zooms the distance between attracting
and repelling cycles for several currents. The reduced-HH gating and
fixed-point/Jacobian helpers shared by the last five live in one cell
(`hh_reduced_*`).

## What to look for

First check how well `ERISIR_REDUCED` follows its full counterpart. For HH,
read `HH_REDUCED_FIXED_POINTS` together with `HH_REDUCED_FP_EVS`; a real-part
crossing signals the stability change. The two cycle examples then show the
nearby attracting and repelling invariant curves.

## Suggested order

1. Run `ERISIR_REDUCED` and `ERISIR_2D_FP`.
2. Run the three HH fixed-point and eigenvalue examples.
3. Finish with `HH_REDUCED_REPELLING_CYCLE` and `HH_REDUCED_CYCLE_DISTANCE`.

## Prerequisites and related chapters

Chapter 10 provides reduced HH geometry and Chapter 13 provides Hopf normal
forms. Chapter 17 compares the finite-onset frequency of these models with
other excitability types.

## Running the examples

Open [`chapter14.ipynb`](chapter14.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom.
