# Bursting

## Overview

These examples create bursts by adding a slow potassium current to Erisir and
INaP-I$_K$ neurons. The slow state moves the fast subsystem between silent
and spiking attractors, yielding alternation between quiescence and rapid
spikes.

## Core ideas

Bursting is a slow-fast rhythm: for an almost fixed slow gate, the fast
subsystem has either a resting state or a spiking cycle. As the slow potassium
gate rises during spiking, it reduces effective drive; as it falls during
silence, it restores excitability. Hysteresis between the fast-subsystem
transitions closes the loop.

## Essential model

The added slow outward current is

$$
I_{K,\rm slow}=g_{K,\rm slow}n_{\rm slow}(v_K-v),\qquad
\dot n_{\rm slow}=\frac{n_{\rm slow,\infty}(v)-n_{\rm slow}}
{\tau_{n,\rm slow}}.
$$

Here $n_{\rm slow}$ is the slow potassium activation, $g_{K,\rm slow}$
is its conductance, $v_K$ is the potassium reversal potential, and
$\tau_{n,\rm slow}$ separates its evolution from the fast voltage and gate
dynamics.

## Code examples

- [`ERISIR_PLUS_SLOW_I_K`](ERISIR_PLUS_SLOW_I_K/) simulates an Erisir burst
  with slow potassium feedback in `fig.png`.
- [`ERISIR_SHOW_SLOW_I_K`](ERISIR_SHOW_SLOW_I_K/) plots Erisir voltage with the
  resulting effective drive and transition levels in `fig.png`.
- [`ELLIPSES`](ELLIPSES/) marks selected portions of the Erisir bursting trace
  in `fig.png`.
- [`INAPIK_PLUS_SLOW_I_K`](INAPIK_PLUS_SLOW_I_K/) simulates INaP-I$_K$
  bursting with the reference slow conductance in `fig.png`.
- [`INAPIK_PLUS_WEAK_SLOW_I_K`](INAPIK_PLUS_WEAK_SLOW_I_K/) shows the response
  with weaker slow potassium feedback in `fig.png`.
- [`INAPIK_PLUS_STRONG_SLOW_I_K`](INAPIK_PLUS_STRONG_SLOW_I_K/) shows the
  response with stronger slow potassium feedback in `fig.png`.
- [`INAPIK_SHOW_SLOW_I_K`](INAPIK_SHOW_SLOW_I_K/) plots INaP-I$_K$ voltage
  and effective drive thresholds in `fig.png`.
- [`SQUARE_WAVES`](SQUARE_WAVES/) marks square-wave burst epochs in `fig.png`.
- [`INAPIK_PLUS_SLOW_I_K_3D`](INAPIK_PLUS_SLOW_I_K_3D/) plots the settled
  INaP-I$_K$ cycle in $(v,n,n_{\rm slow})$ space in `fig.png`.

## What to look for

Compare the weak, reference, and strong INaP-I$_K$ slow conductances to see
how slow feedback changes the rhythm. The two `SHOW_SLOW_I_K` plots make the
effective current cross thresholds while the voltage alternates. Use
`INAPIK_PLUS_SLOW_I_K_3D` to connect the trace to a three-dimensional loop.

## Suggested order

1. Run `INAPIK_PLUS_SLOW_I_K`, its weak and strong variants, and `SQUARE_WAVES`.
2. Run `INAPIK_SHOW_SLOW_I_K` and `INAPIK_PLUS_SLOW_I_K_3D`.
3. Compare `ERISIR_PLUS_SLOW_I_K`, `ERISIR_SHOW_SLOW_I_K`, and `ELLIPSES`.

## Prerequisites and related chapters

Chapter 15 introduces slow-fast canard behavior, Chapter 17 provides the
fast-subsystem f--I interpretation, and Chapter 18 shows how slow currents
create hysteresis and bistability.

## Running the examples

Run `python main.py` from each immediate example directory. The scripts use
NumPy and Matplotlib and save `fig.png`; the long Erisir and 3D trajectories
may take longer than the short INaP-I$_K$ simulations.
