# Bursting

## Overview

These examples create bursts by adding a slow potassium current to Erisir
and INaP-I$_K$ neurons. The slow state moves the fast subsystem between
silent and spiking attractors, yielding alternation between quiescence and
rapid spikes.

## Core ideas

Bursting is a slow-fast rhythm: for an almost fixed slow gate, the fast
subsystem has either a resting state or a spiking cycle. As the slow
potassium gate rises during spiking, it reduces effective drive; as it
falls during silence, it restores excitability. Hysteresis between the
fast-subsystem transitions closes the loop.

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

All nine examples now live in one notebook, [`chapter19.ipynb`](chapter19.ipynb):
`simulate_inapik_plus_slow_i_k` simulates an INaP-I$_K$ burst with the
reference slow conductance ($g_{K,\rm slow}=5$; an interactive slider lets
you sweep it); `simulate_inapik_plus_weak_slow_i_k`/
`simulate_inapik_plus_strong_slow_i_k` show the response with weaker/stronger
slow potassium feedback; `simulate_inapik_show_slow_i_k` plots INaP-I$_K$
voltage with the effective drive and transition thresholds;
`simulate_square_waves` marks the quasi-steady (square-wave) portions of an
INaP-I$_K$ burst; `simulate_inapik_plus_slow_i_k_3d` plots the settled
INaP-I$_K$ cycle in $(v,n,n_{\rm slow})$ space; `simulate_erisir_plus_slow_i_k`
simulates an Erisir burst with slow potassium feedback and is reused by
`plot_erisir_show_slow_i_k` (voltage with effective drive and thresholds)
and `plot_ellipses` (local extrema marked on the trace).

## What to look for

Compare the weak, reference, and strong INaP-I$_K$ slow conductances (or
sweep `g_k_slow` with the interactive slider) to see how slow feedback
changes the rhythm. The `_show_slow_i_k` plots make the effective current
cross thresholds while the voltage alternates. Use
`simulate_inapik_plus_slow_i_k_3d` to connect the trace to a
three-dimensional loop.

## Suggested order

1. Run `simulate_inapik_plus_slow_i_k`, its weak and strong variants, and
   `simulate_square_waves`.
2. Run `simulate_inapik_show_slow_i_k` and `simulate_inapik_plus_slow_i_k_3d`.
3. Compare `simulate_erisir_plus_slow_i_k` with `plot_erisir_show_slow_i_k`
   and `plot_ellipses`.

## Prerequisites and related chapters

Chapter 15 introduces slow-fast canard behavior, Chapter 17 provides the
fast-subsystem f--I interpretation, and Chapter 18 shows how slow currents
create hysteresis and bistability.

## Running the examples

Open [`chapter19.ipynb`](chapter19.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom.
