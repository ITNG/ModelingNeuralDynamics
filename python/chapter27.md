# Phase locking with delays

## Overview

Propagation delays change the phase at which an oscillator receives a pulse,
so they change the locking map itself. These examples simulate delayed pulse
coupling for two and three oscillators and realize delayed locking with theta
neurons.

## Core ideas

A delayed pulse arrives after its source phase has advanced, and this timing
must be included before applying a PRC. Consequently a delay can stabilize a
phase difference that was not locked without delay, or destabilize synchronous
timing. Three oscillators add consistent pairwise timing constraints; their
locked pattern need not be a simple two-cell extension.

## Essential model

For natural period $T$ and transmission delay $d$, an event map applies a
phase shift at a phase displaced by $d/T$. In a theta-neuron realization, the
continuous angle evolves between discrete arrivals, so the delay is held as
an event-time condition rather than being folded into a static voltage term.

## Code examples

All three examples now live in one notebook, [`chapter27.ipynb`](chapter27.ipynb):
`simulate_two_delayed_pulse_coupled_osc` iterates a two-oscillator phase map
with delayed pulses, running a short-delay ($\delta=0.1$) and a long-delay
($\delta=0.7$) case with the shared `simulate_two_delayed_pulse_pair` engine
(`two_pulse_g`/`two_pulse_f` give the pulse and interaction-function maps).
`simulate_three_delayed_pulse_coupled_osc` extends delayed pulse coupling to
three all-to-all oscillators using `simulate_three_delayed_pulse_pair`
(`three_pulse_g` for the pulse map), comparing $\delta=0.45$ against
$\delta=0.55$ and returning each run's event times for
`plot_three_delayed_pulse_coupled_osc`. `simulate_two_theta_neurons_grid`
classifies a $9\times9$ grid of delayed, pulse-coupled theta-neuron
simulations (`simulate_two_theta_neurons_pair`, `theta_neuron_inc`,
`sync_measure`) as synchronized or unsynchronized and
`plot_two_theta_neurons_grid` draws that region in the $(\epsilon,\delta)$
parameter plane.

## What to look for

Vary the initial conditions mentally while following the event sequence: a
locked state repeats the same pulse-arrival phases. In the three-cell plot,
distinguish a repeating collective order from exact simultaneous spikes. For
the theta pair, read the red and blue grid points as synchronized and
unsynchronized outcomes, respectively, and compare them with the plotted
boundary in the $(\epsilon,\delta)$ plane.

## Suggested order

1. Run `simulate_two_delayed_pulse_coupled_osc`.
2. Run `simulate_three_delayed_pulse_coupled_osc` and identify its repeating
   order.
3. Run `simulate_two_theta_neurons_grid` to connect the phase-map result to
   continuous neuron dynamics.

## Prerequisites and related chapters

Chapter 26 establishes undelayed pulse-coupled phase maps, and Chapter 8
introduces theta neurons. Chapter 28 treats weak coupling, while Chapter 29
focuses on the stability of a synchronous phase relation.

## Running the examples

Open [`chapter27.ipynb`](chapter27.ipynb) in Jupyter, or via the Colab badge
at the top of the notebook, and run all cells top to bottom. The theta-neuron
grid cell is the slowest (about a minute and a half): it runs 81 sampled
$(\epsilon,\delta)$ pairs, each a long fine-step theta-neuron simulation.
