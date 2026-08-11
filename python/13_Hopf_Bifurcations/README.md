# Hopf bifurcations

## Overview

These normal-form examples contrast supercritical and subcritical Hopf
bifurcations through radial equations, phase planes, and bifurcation diagrams.

## Core ideas

A Hopf bifurcation changes the stability of an equilibrium as a complex pair
of eigenvalues crosses the imaginary axis. A supercritical Hopf creates a
small stable cycle, whereas a subcritical Hopf has an unstable small cycle and
can coexist with a larger attracting cycle. Solid and dashed branches in the
plots distinguish attracting and repelling structures.

## Essential model

With polar radius $r$, the normal forms used here are

$$
\dot r=Ir-r^3 \quad\text{(supercritical)},\qquad
\dot r=Ir+r^3 \quad\text{(subcritical)}.
$$

Here $I$ is the bifurcation parameter and $r\ge0$ is distance from the
equilibrium. The extended subcritical examples add a $-r^5$ term, which
supports both an unstable inner and a stable outer cycle.

## Code examples

- [`HOPF_SUP`](HOPF_SUP/) plots the supercritical radial vector field in
  `fig.png`.
- [`HOPF_SUP_BIF_DIAG`](HOPF_SUP_BIF_DIAG/) draws its equilibrium and stable
  cycle branches in `fig.png`.
- [`HOPF_SUP_PHASE_PLANE`](HOPF_SUP_PHASE_PLANE/) integrates supercritical
  spirals and cycles in `fig.png`.
- [`HOPF_SUB`](HOPF_SUB/) plots the subcritical radial vector field in
  `fig.png`.
- [`HOPF_SUB_BIF_DIAG`](HOPF_SUB_BIF_DIAG/) draws the subcritical unstable
  cycle branch in `fig.png`.
- [`HOPF_SUB_PHASE_PLANE`](HOPF_SUB_PHASE_PLANE/) shows the corresponding
  repelling cycle and phase-plane trajectories in `fig.png`.
- [`HOPF_SUB_2`](HOPF_SUB_2/) adds the quintic saturation term in `fig.png`.
- [`HOPF_SUB_BIF_DIAG_2`](HOPF_SUB_BIF_DIAG_2/) plots its inner repelling and
  outer attracting cycles in `fig.png`.
- [`HOPF_SUB_PHASE_PLANE_2`](HOPF_SUB_PHASE_PLANE_2/) traces the two-cycle
  phase-plane organization in `fig.png`.

## What to look for

Compare `HOPF_SUP_BIF_DIAG` with `HOPF_SUB_BIF_DIAG`: the stable cycle appears
on opposite sides of onset. Then use the phase-plane scripts to connect a
branch's line style to trajectories moving toward or away from its cycle.

## Suggested order

1. Run the three `HOPF_SUP` examples.
2. Run `HOPF_SUB`, `HOPF_SUB_BIF_DIAG`, and `HOPF_SUB_PHASE_PLANE`.
3. Use the three `_2` examples to study coexistence of two cycles.

## Prerequisites and related chapters

Chapter 11 introduces local bifurcations and Chapter 12 uses Jacobian
eigenvalues. Chapters 14--15 apply Hopf and subcritical-cycle geometry to
conductance-based neuron reductions.

## Running the examples

Run `python main.py` from any immediate example directory. The scripts use
NumPy and Matplotlib and save `fig.png` locally.
