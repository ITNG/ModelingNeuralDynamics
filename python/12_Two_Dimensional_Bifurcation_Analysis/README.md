# Two-dimensional bifurcation analysis

## Overview

The historical folder title is retained, but these scripts correspond to the
book's type-1 onset material. They reduce the RTM neuron to two variables and
show its fixed points, invariant cycles, and saddle-node-on-invariant-circle
(SNIC) onset geometry.

## Core ideas

In a two-dimensional reduction, fixed points are intersections of nullclines
and their Jacobian eigenvalues classify stability. A periodic orbit is an
invariant cycle. At type-1 onset, a saddle and node meet on the cycle, so the
orbit spends an increasingly long time near the collision and its firing rate
approaches zero.

## Essential model

The reduced RTM dynamics use instantaneous activation and the approximation
$m=m_\infty(v)$, $h=1-n$:

$$
C\dot v=I_{\rm Na}(v,n)+I_{\rm K}(v,n)+I_L(v)+I,\qquad
\dot n=\alpha_n(v)(1-n)-\beta_n(v)n.
$$

Here $v$ is membrane voltage, $n$ is potassium activation, $I$ is
applied current, and the sodium, potassium, and leak terms use RTM
conductances and reversal potentials.

## Code examples

- [`RTM_2D_FP`](RTM_2D_FP/) scans applied current, finds reduced RTM fixed
  points, classifies their Jacobian eigenvalues, and saves `fig.png`.
- [`RTM_2D_INVARIANT_CYCLE`](RTM_2D_INVARIANT_CYCLE/) plots trajectories and
  crossings around the RTM invariant cycle and its onset geometry in `fig.png`.

## What to look for

`RTM_2D_FP` separates branches by stability and eigenvalue type. In
`RTM_2D_INVARIANT_CYCLE`, compare trajectories at and around the critical
current: motion slows near the saddle-node region, which is the phase-plane
signature behind the low-frequency type-1 onset.

## Suggested order

1. Run `RTM_2D_FP` to locate and classify equilibria.
2. Run `RTM_2D_INVARIANT_CYCLE` to see how trajectories use those structures.
3. Continue to Chapter 17 to measure the resulting frequency-current curve.

## Prerequisites and related chapters

Chapter 10 introduces reductions and nullclines, and Chapter 11 introduces
saddle-node collisions. Chapter 14 contrasts type-2 onset; Chapter 17 turns
these onset mechanisms into firing-rate curves.

## Running the examples

Run `python main.py` from either immediate example directory. Each script uses
NumPy and Matplotlib, and `RTM_2D_INVARIANT_CYCLE` also uses the local arrow
plotting helper; both save `fig.png`.
