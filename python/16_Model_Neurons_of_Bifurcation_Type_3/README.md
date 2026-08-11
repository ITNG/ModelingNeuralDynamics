# Model neurons of bifurcation type 3

## Overview

This chapter presents type-3, or phasic, excitability. Persistent sodium and
potassium dynamics generate transient responses, and self-exciting theta
neurons show the same geometry in a phase model with a slow feedback state.

## Core ideas

A type-3 neuron can answer a step with only a brief spike or a few spikes,
rather than tonic firing. Fixed points and phase portraits organize this
phasic response. In a self-exciting theta neuron, a slow state jumps or rises
around spikes and changes the phase velocity, creating a geometric feedback
mechanism.

## Essential model

The persistent-sodium/potassium (INaP-I$_K$) reduction uses

$$
C\dot v=g_{\rm Na}m_\infty(v)(v_{\rm Na}-v)+g_Kn(v_K-v)
+g_L(v_L-v)+I,\qquad \dot n=(n_\infty(v)-n)/\tau_n.
$$

Here $v$ is voltage, $n$ is potassium activation, $m_\infty$ is
instantaneous persistent sodium activation, and $I$ is applied current.

## Code examples

- [`INAPIK_FIXED_POINTS`](INAPIK_FIXED_POINTS/) scans current and classifies
  INaP-I$_K$ equilibria in `fig.png`.
- [`INAPIK_PHASE_PLANE`](INAPIK_PHASE_PLANE/) plots phase-plane trajectories,
  fixed points, and cycles at several currents in `fig.png`.
- [`SELF_EXCITING_THETA_NEURON`](SELF_EXCITING_THETA_NEURON/) simulates a
  theta neuron with discrete slow-state increments and saves `fig.png`.
- [`SELF_EXCITING_THETA_SMOOTH`](SELF_EXCITING_THETA_SMOOTH/) uses smooth slow
  feedback in the theta model and saves `fig.png`.
- [`SETN_PHASE_PLANE`](SETN_PHASE_PLANE/) draws theta--slow-state phase planes
  and threshold markers in `fig.png`.

## What to look for

Use `INAPIK_FIXED_POINTS` before the multi-panel `INAPIK_PHASE_PLANE` so the
markers have a clear meaning. Then compare the jump-based and smooth theta
feedback in the two `SELF_EXCITING` examples; `SETN_PHASE_PLANE` reveals the
threshold structure underlying those traces.

## Suggested order

1. Run `INAPIK_FIXED_POINTS` and `INAPIK_PHASE_PLANE`.
2. Run `SELF_EXCITING_THETA_NEURON` and its `_SMOOTH` variant.
3. Use `SETN_PHASE_PLANE` to connect the time traces to phase geometry.

## Prerequisites and related chapters

Chapter 08 introduces theta neurons. Chapters 10--14 supply phase-plane and
bifurcation language, and Chapter 17 compares the firing-rate response of the
INaP-I$_K$ and theta models.

## Running the examples

Run `python main.py` from each immediate example directory. Every script uses
NumPy and Matplotlib and saves `fig.png` locally.
