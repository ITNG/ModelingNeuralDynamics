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

Let \(P_0\) and \(P_1\) denote the two pulse-to-spike timing maps used to
construct a synchronous-state return map \(P\). A synchronous timing relation
is locally contracting when the magnitude of its derivative is below one. The
condition-number calculations measure how perturbations to the maps affect this
conclusion.

## Code examples

- [`ILLUSTRATE_P0_AND_P1`](ILLUSTRATE_P0_AND_P1/) traces trajectories used to
  illustrate the two timing maps before and after a pulse.
- [`LIF_P_AND_S`](LIF_P_AND_S/) computes and plots the LIF return map together
  with its stability-related quantities.
- [`LIF_CONDITION_NUMBERS`](LIF_CONDITION_NUMBERS/) evaluates LIF condition
  numbers for the pulse-map calculation.
- [`LIF_WITH_INHIBITORY_PULSE`](LIF_WITH_INHIBITORY_PULSE/) simulates LIF cells
  receiving a common inhibitory pulse and records their timing response.
- [`RTM_WITH_INHIBITORY_PULSE`](RTM_WITH_INHIBITORY_PULSE/) performs the
  corresponding conductance-based RTM inhibitory-pulse experiment.
- [`RTM_CONDITION_NUMBERS`](RTM_CONDITION_NUMBERS/) evaluates sensitivity of
  the RTM synchronous-state calculation.
- [`RIVER`](RIVER/) plots the theta-neuron river geometry that organizes how
  nearby trajectories are drawn toward or away from a timing relation.

## What to look for

In the pulse simulations, compare the timing gap before and after the common
inhibition rather than just the voltage deflection. Read the LIF and RTM
condition-number plots alongside their maps: a contraction conclusion with a
large condition number deserves more caution. In `RIVER`, follow neighbouring
theta trajectories to connect the geometric flow with phase contraction.

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
