# Chemical synapses

## Overview

Chemical synapses turn presynaptic activity into a conductance that drives the
postsynaptic voltage toward a reversal potential. These examples use RTM
neurons to make the synaptic gate visible, then use that gate for a self-synapse
and for the temporal accumulation produced by repeated input.

## Core ideas

The synaptic current is conductance based:
$I_{\rm syn}=g_{\rm syn}s(v_{\rm syn}-v)$. A release variable can rise rapidly
after a spike and feed
the gate $s$, which decays more slowly. The rise and decay constants therefore
shape both the peak time and the duration of the conductance. The NMDA factor
also depends on postsynaptic voltage because magnesium block is relieved by
depolarization.

## Essential model

For the two-stage synapse used here, release $q$ and activation $s$ evolve
as

$$
\dot q=R(v)(1-q)-q/\tau_{d,q},\qquad
\dot s=q(1-s)/\tau_r-s/\tau_d.
$$

`tau_d_q_function` numerically chooses $\tau_{d,q}$ to obtain a requested
activation peak time. The resulting $s$ multiplies the synaptic conductance
in the RTM voltage equation.

## Code examples

All eight examples now live in one notebook, [`chapter20.ipynb`](chapter20.ipynb):
`simulate_b_jahr_stevens` plots the voltage-dependent NMDA magnesium-block
factor from Jahr and Stevens; `simulate_rtm_plot_s` compares a fast and a
slower synaptic-gate rise alongside the RTM voltage trace;
`simulate_rtm_plot_q` separates transmitter-release $q$ from the gate $s$,
showing why a two-variable synapse has a delayed profile;
`simulate_two_stage_synapse` (shared by three examples below) uses the
explicit release-and-gate system: called directly with fixed time constants
for two timing choices, or with `tau_d_q_function`-solved time constants to
hit a prescribed peak time; `simulate_rtm_with_autapse_f_i_curve` follows
forward and backward RTM frequency-current sweeps with an excitatory
autapse; `simulate_s_buildup` shows how closely spaced presynaptic events
build up a synaptic gate, called once with a fast decay and once with a
slower one (buildup and slow-buildup).

## What to look for

In the gate plots, locate the delay between a voltage spike and the maximum of
$s$. Increasing the rise or release time moves that maximum and broadens the
conductance. Compare the two buildup calls to see that a slow decay retains
activation between events. In the autapse sweep, compare the two directions to
identify any history dependence of the firing state.

## Suggested order

1. Run `simulate_rtm_plot_s`, `simulate_rtm_plot_q`, and
   `simulate_two_stage_synapse` with the two fixed-timing calls.
2. Run the prescribed-peak-time `simulate_two_stage_synapse` calls, then
   `simulate_s_buildup` with the fast and slow decay parameters.
3. Examine `simulate_b_jahr_stevens` and `simulate_rtm_with_autapse_f_i_curve`.

## Prerequisites and related chapters

The RTM conductance model is introduced earlier in the single-neuron chapters.
Chapter 21 replaces chemical conductance with electrical coupling, while
Chapters 23--29 use pulses and phase descriptions to study network timing.

## Running the examples

Open [`chapter20.ipynb`](chapter20.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom. The
`simulate_rtm_with_autapse_f_i_curve` cell integrates a forward+backward
sweep over 31 values of $I$; its inner loop is JIT-compiled with numba, so
after the first (one-time compile) call it takes well under a second
instead of several minutes per direction.
