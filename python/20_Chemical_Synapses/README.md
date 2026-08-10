# Chemical synapses

## Overview

Chemical synapses turn presynaptic activity into a conductance that drives the
postsynaptic voltage toward a reversal potential. These examples use RTM
neurons to make the synaptic gate visible, then use that gate for a self-synapse
and for the temporal accumulation produced by repeated input.

## Core ideas

The synaptic current is conductance based: \(I_{\rm syn}=g_{\rm syn}s
(v_{\rm syn}-v)\). A release variable can rise rapidly after a spike and feed
the gate \(s\), which decays more slowly. The rise and decay constants therefore
shape both the peak time and the duration of the conductance. The NMDA factor
also depends on postsynaptic voltage because magnesium block is relieved by
depolarization.

## Essential model

For the two-stage synapse used here, release \(q\) and activation \(s\) evolve
as

\[
\dot q=R(v)(1-q)-q/\tau_{d,q},\qquad
\dot s=q(1-s)/\tau_r-s/\tau_d.
\]

`tau_d_q_function` numerically chooses \(\tau_{d,q}\) to obtain a requested
activation peak time. The resulting \(s\) multiplies the synaptic conductance
in the RTM voltage equation.

## Code examples

- [`B_JAHR_STEVENS`](B_JAHR_STEVENS/) plots the voltage-dependent NMDA
  magnesium-block factor from Jahr and Stevens.
- [`RTM_PLOT_S`](RTM_PLOT_S/) compares a fast and a slower synaptic-gate rise
  alongside the RTM voltage trace.
- [`RTM_PLOT_Q`](RTM_PLOT_Q/) separates transmitter-release \(q\) from the
  gate \(s\), showing why a two-variable synapse has a delayed profile.
- [`RTM_PLOT_S_TWO_VARIABLES`](RTM_PLOT_S_TWO_VARIABLES/) uses the explicit
  release-and-gate system for two timing choices.
- [`RTM_PLOT_S_PRESCRIBE_TAU_PEAK`](RTM_PLOT_S_PRESCRIBE_TAU_PEAK/) solves for
  the release decay constant that gives each prescribed peak time.
- [`RTM_WITH_AUTAPSE_F_I_CURVE`](RTM_WITH_AUTAPSE_F_I_CURVE/) follows forward
  and backward RTM frequency--current sweeps with an excitatory autapse.
- [`S_BUILDUP`](S_BUILDUP/) shows how closely spaced presynaptic events build
  up a synaptic gate.
- [`S_SLOW_BUILDUP`](S_SLOW_BUILDUP/) repeats the buildup experiment with a
  slower synaptic time course.

## What to look for

In the gate plots, locate the delay between a voltage spike and the maximum of
\(s\). Increasing the rise or release time moves that maximum and broadens the
conductance. Compare the two buildup plots to see that a slow decay retains
activation between events. In the autapse sweep, compare the two directions to
identify any history dependence of the firing state.

## Suggested order

1. Run `RTM_PLOT_S`, `RTM_PLOT_Q`, and `RTM_PLOT_S_TWO_VARIABLES`.
2. Run `RTM_PLOT_S_PRESCRIBE_TAU_PEAK`, then `S_BUILDUP` and `S_SLOW_BUILDUP`.
3. Examine `B_JAHR_STEVENS` and `RTM_WITH_AUTAPSE_F_I_CURVE`.

## Prerequisites and related chapters

The RTM conductance model is introduced earlier in the single-neuron chapters.
Chapter 21 replaces chemical conductance with electrical coupling, while
Chapters 23--29 use pulses and phase descriptions to study network timing.

## Running the examples

Run `python main.py` from an immediate example directory. The scripts use
NumPy, SciPy, and Matplotlib and save their figures locally; the autapse sweep
can take longer because it integrates each input level to steady state.
