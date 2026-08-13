# Stability of the synchronous state

## Overview

Synchrony is stable when an inhibitory perturbation contracts a small timing
difference rather than amplifying it. These examples calculate the
contraction maps and their condition numbers for LIF and RTM cells, then
give a geometric theta-neuron "river" picture of the same sensitivity.

## Core ideas

Two nearby phases can be represented by their pre-pulse and post-pulse time
separations. An inhibitory pulse may contract that separation, but numerical
or parameter sensitivity can make the inferred stability fragile. Condition
numbers quantify this sensitivity. LIF reset dynamics and smooth RTM
dynamics can therefore react differently even when both receive an
inhibitory pulse.

## Essential model

Let $P_0$ and $P_1$ denote the two pulse-to-spike timings. The LIF examples
define their mean timing and normalized separation as

$$
P=\frac{P_0+P_1}{2},\qquad S=\frac{P_0-P_1}{P}.
$$

Here $S$, rather than a derivative of $P$, is the plotted stability
quantity. The sensitivity calculations report how changes in parameters
alter the mean timing and its inference.

## Code examples

All seven examples now live in one notebook, [`chapter29.ipynb`](chapter29.ipynb):
`simulate_p0_and_p1` traces the two timing maps before and after a pulse
(`plot_p0_and_p1`). `P0`, `P1`, `P`, and `S` are the shared LIF timing maps
used by both LIF examples below: `simulate_lif_p_and_s` plots the LIF mean
timing $P$ and normalized separation $S$ while varying synaptic decay,
conductance, and drive (`plot_lif_p_and_s`); `compute_lif_condition_numbers`
returns a dictionary of baseline LIF mean timings and percent changes under
parameter perturbations for a weak/slow and a strong/fast synapse.
`simulate_lif_pulse_panels` simulates LIF cells receiving a common
inhibitory pulse and records their timing response (`plot_lif_pulse_panels`).
`rtm_init` and `simulate_rtm` are shared by the two RTM examples:
`simulate_rtm_pulse_panels` performs the RTM analogue of the LIF
inhibitory-pulse experiment (`plot_rtm_pulse_panels`, reused by the
condition-numbers example below); `compute_rtm_condition_numbers` returns
RTM timing sensitivities (and the same three voltage traces) for no,
weak/slow, and strong/fast inhibitory conductances. `simulate_river` and
`plot_river` draw the theta-neuron river geometry that organizes how nearby
trajectories are drawn toward or away from a timing relation.

## What to look for

In the pulse simulations, compare the timing gap before and after the
common inhibition rather than just the voltage deflection. Read the printed
LIF and RTM sensitivity results alongside the $P$ and $S$ plots and the RTM
voltage traces: a conclusion that changes strongly with parameters deserves
more caution. In the river plot, follow neighbouring theta trajectories to
connect the geometric flow with phase contraction.

## Suggested order

1. Run `simulate_p0_and_p1`, `simulate_lif_p_and_s`, and
   `simulate_lif_pulse_panels`.
2. Run `compute_lif_condition_numbers` and compare the sensitivity with the
   LIF map.
3. Run `simulate_rtm_pulse_panels`, `compute_rtm_condition_numbers`, and
   `simulate_river`.

## Prerequisites and related chapters

Chapter 7 provides LIF reset dynamics, Chapter 25 supplies PRCs, and
Chapters 26-28 introduce phase-map, delay, and weak-coupling stability
ideas. The theta-neuron geometry builds on Chapter 8.

## Running the examples

Open [`chapter29.ipynb`](chapter29.ipynb) in Jupyter, or via the Colab badge
at the top of the notebook, and run all cells top to bottom. RTM
calculations and condition-number sweeps take a little longer than the LIF
and geometric examples but still finish in seconds.
