# Phase response curves (PRCs)

## Overview

A phase response curve records how a perturbation applied at one point of an
oscillation advances or delays the next cycle. These examples compute finite
and weak-pulse PRCs for RTM and WB neurons, compare them with theta-neuron
formulae, and turn a PRC into an interaction function.

## Core ideas

A finite PRC measures a measurable phase change from a pulse; an infinitesimal
PRC is its small-amplitude limit. A type-1 PRC is predominantly advancing,
whereas a type-2 PRC has advancing and delaying regions. Pulse duration and
strength matter: a brief voltage kick, a finite synaptic pulse, and a weaker
perturbation need not produce the same curve. Averaging the response over a
cycle yields an interaction function for weakly coupled oscillators.

## Essential model

For unperturbed period \(T\), let a pulse applied at phase \(\phi\) change the
next-cycle timing by \(\Delta t\). A phase response can be represented as

\[
Z(\phi)=-\Delta t/T.
\]

For weak coupling, the interaction function combines a phase response and a
synaptic waveform; it describes how one oscillator's phase affects another's
slow phase drift.

## Code examples

- [`MISC_PRC`](MISC_PRC/) collects additional finite-pulse PRC constructions
  and comparison plots.
- [`PHASE_SHIFT`](PHASE_SHIFT/) measures the RTM time shift after a pulse at a
  selected phase.
- [`RTM_PRC`](RTM_PRC/) samples a finite synaptic-pulse PRC for RTM neurons.
- [`RTM_PRC_SHORT`](RTM_PRC_SHORT/) uses an instantaneous short voltage kick
  to obtain an RTM PRC.
- [`RTM_PRC_WEAK`](RTM_PRC_WEAK/) evaluates the weak-pulse RTM response.
- [`RTM_PRC_SHORT_AND_WEAK`](RTM_PRC_SHORT_AND_WEAK/) varies small kick sizes
  to expose the infinitesimal-response limit.
- [`RTM_PRC_THREE_WEAK_ONES`](RTM_PRC_THREE_WEAK_ONES/) compares three weak
  RTM synaptic perturbations.
- [`RTM_INTERACTION_FUNCTION`](RTM_INTERACTION_FUNCTION/) derives and plots an
  RTM interaction function from pulse responses.
- [`THETA_PRC`](THETA_PRC/) plots the theta-neuron PRC.
- [`THETA_PRC_SHORT_WEAK`](THETA_PRC_SHORT_WEAK/) shows the theta response to
  short, weak inputs on a dense phase grid.
- [`THETA_F`](THETA_F/) plots the theta-neuron phase map induced by a pulse.
- [`THETA_F_TILDE`](THETA_F_TILDE/) plots the associated transformed phase map.
- [`WB_PRC_INHIBITORY_PULSE`](WB_PRC_INHIBITORY_PULSE/) measures a WB PRC for
  an inhibitory synaptic pulse.

## What to look for

Read positive and negative portions as advances and delays using the plotted
sign convention. Compare `RTM_PRC` with its short and weak variants before
calling a curve infinitesimal. The theta plots provide a compact analytic
reference; the WB inhibitory curve makes clear that pulse polarity changes the
timing effect. In `RTM_INTERACTION_FUNCTION`, identify zeros because they are
candidate phase relationships for coupled cells.

## Suggested order

1. Run `PHASE_SHIFT`, `RTM_PRC`, and `RTM_PRC_SHORT`.
2. Compare `RTM_PRC_WEAK`, `RTM_PRC_SHORT_AND_WEAK`, and
   `RTM_PRC_THREE_WEAK_ONES`.
3. Run the `THETA_*` examples, `WB_PRC_INHIBITORY_PULSE`, and
   `RTM_INTERACTION_FUNCTION`; use `MISC_PRC` for further comparisons.

## Prerequisites and related chapters

Chapter 8 introduces theta neurons, and Chapter 23 motivates phase maps from
pulsed input. Chapters 26--28 use PRCs and interaction functions to determine
locking and phase differences in coupled oscillators.

## Running the examples

Run `python main.py` in an immediate example directory. The numerical PRC
scripts use NumPy and Matplotlib and may integrate many initial phases, so they
are slower than the theta formula plots.
