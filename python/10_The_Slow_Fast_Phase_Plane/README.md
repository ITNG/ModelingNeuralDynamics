# The slow-fast phase plane

## Overview

These examples reduce conductance-based spiking to phase-plane geometry. They
compare the visibly separated slow and fast variables in FitzHugh--Nagumo and
Hodgkin--Huxley (HH) systems, then show how that separation organizes a
periodic orbit.

## Core ideas

When one state changes much faster than another, trajectories move quickly
between branches of a slow manifold and linger near slow branches. Nullclines
locate the directions with zero velocity, while a closed orbit describes
repeated firing. A reduction with instantaneous sodium activation and
$h+n=0.83$ makes this geometry visible in a two-dimensional HH system.

## Essential model

The FitzHugh--Nagumo example uses

$$
\dot v=v-v^3/3-n+I,\qquad \dot n=(av-n)/\tau_n.
$$

Here $v$ is the fast voltage-like variable, $n$ is the slow recovery
variable, $I$ is applied current, $a$ sets the recovery nullcline, and
$\tau_n$ makes recovery slow. The reduced HH examples take
$m=m_\infty(v)$ and $h=0.83-n$, leaving $(v,n)$ as the phase plane.

## Code examples

- [`FN`](FN/) integrates FitzHugh--Nagumo and plots its nullclines, orbit, and
  voltage trace in `fig.png`.
- [`HH_H_PLUS_N`](HH_H_PLUS_N/) plots the HH combination $h+n$ against its
  `0.83` approximation in `fig_10_1.png`.
- [`REDUCED_HH`](REDUCED_HH/) compares the reduced HH voltage and gates with
  the imposed $h=0.83-n$ relation in `fig_10_2.png`.
- [`HH_NULLCLINES_PLUS_SOLUTION`](HH_NULLCLINES_PLUS_SOLUTION/) draws the
  reduced HH nullclines and a trajectory with direction arrows in `fig.png`.
- [`HH_CYCLE_SPEED`](HH_CYCLE_SPEED/) colors portions of the HH orbit by their
  phase-plane speed in `fig.png`.

## What to look for

In `FN`, find the long motion near the cubic nullcline and the rapid jumps
between its outer branches. Compare `HH_H_PLUS_N` with `REDUCED_HH` before
using the two-dimensional HH plots. In `HH_CYCLE_SPEED`, green points mark
the slow portions of the orbit and blue points the faster traversal.

## Suggested order

1. Run `FN` to connect nullclines to slow-fast motion.
2. Run `HH_H_PLUS_N` and `REDUCED_HH` to inspect the HH reduction.
3. Use `HH_NULLCLINES_PLUS_SOLUTION` and `HH_CYCLE_SPEED` to read its cycle.

## Prerequisites and related chapters

This chapter builds on HH gates from Chapters 03--04 and the reduced models
introduced in Chapter 05. Chapters 11--15 use the same phase-plane language
to study bifurcations and canards.

## Running the examples

From any immediate example directory, run `python main.py`. The scripts use
NumPy, SciPy where noted, and Matplotlib, and save their named PNG files in
the working directory.
