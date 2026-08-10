# Stability of the synchronous state

## Overview

Synchrony is stable when an inhibitory perturbation contracts a small timing
difference rather than amplifying it. These examples calculate the contraction
maps and their condition numbers for LIF and RTM cells, then give a geometric
theta-neuron ``river'' picture of the same sensitivity.

## Core ideas

Two nearby phases can be represented by their pre-pulse and post-pulse time
separations. An inhibitory pulse may contract that separation, but numerical
or parameter sensitivity can make the inferred stability fragile. Condition
numbers quantify this sensitivity. LIF reset dynamics and smooth RTM dynamics
can therefore react differently even when both receive an inhibitory pulse.

## Essential model

Let \(P_0\) and \(P_1\) denote the two pulse-to-spike timings. The LIF scripts
define their mean timing and normalized separation as

\[
P=\frac{P_0+P_1}{2},\qquad S=\frac{P_0-P_1}{P}.
\]

Here \(S\), rather than a derivative of \(P\), is the plotted stability
quantity. The sensitivity calculations report how changes in parameters alter
the mean timing and its inference.

## Code examples

- [`ILLUSTRATE_P0_AND_P1`](ILLUSTRATE_P0_AND_P1/) traces trajectories used to
  illustrate the two timing maps before and after a pulse.
- [`LIF_P_AND_S`](LIF_P_AND_S/) plots the LIF mean timing \(P\) and normalized
  separation \(S\) while varying synaptic decay, conductance, and drive.
- [`LIF_CONDITION_NUMBERS`](LIF_CONDITION_NUMBERS/) prints a dictionary of
  baseline LIF mean timings and percent changes under parameter perturbations.
- [`LIF_WITH_INHIBITORY_PULSE`](LIF_WITH_INHIBITORY_PULSE/) simulates LIF cells
  receiving a common inhibitory pulse and records their timing response.
- [`RTM_WITH_INHIBITORY_PULSE`](RTM_WITH_INHIBITORY_PULSE/) performs the
  corresponding conductance-based RTM inhibitory-pulse experiment.
- [`RTM_CONDITION_NUMBERS`](RTM_CONDITION_NUMBERS/) prints RTM timing
  sensitivities and plots voltage traces for no, weak/slow, and strong/fast
  inhibitory conductances.
- [`RIVER`](RIVER/) plots the theta-neuron river geometry that organizes how
  nearby trajectories are drawn toward or away from a timing relation.

## What to look for

In the pulse simulations, compare the timing gap before and after the common
inhibition rather than just the voltage deflection. Read the printed LIF and
RTM sensitivity results alongside the \(P\) and \(S\) plots and the RTM voltage
traces: a conclusion that changes strongly with parameters deserves more
caution. In `RIVER`, follow neighbouring theta trajectories to connect the
geometric flow with phase contraction.

## Suggested order

1. Run `ILLUSTRATE_P0_AND_P1`, `LIF_P_AND_S`, and
   `LIF_WITH_INHIBITORY_PULSE`.
2. Run `LIF_CONDITION_NUMBERS` and compare the sensitivity with the LIF map.
3. Run `RTM_WITH_INHIBITORY_PULSE`, `RTM_CONDITION_NUMBERS`, and `RIVER`.

## Prerequisites and related chapters

Chapter 7 provides LIF reset dynamics, Chapter 25 supplies PRCs, and
Chapters 26--28 introduce phase-map, delay, and weak-coupling stability ideas.
The theta-neuron geometry builds on Chapter 8.

## Running the examples

Run `python main.py` in an immediate example directory. The scripts use NumPy
and Matplotlib; RTM calculations and condition-number sweeps can take longer
than the LIF and geometric examples.
