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
phase shift at a phase displaced by $d/T$. In a theta-neuron realization,
the continuous angle evolves between discrete arrivals, so the delay is held as
an event-time condition rather than being folded into a static voltage term.

## Code examples

- [`TWO_DELAYED_PULSE_COUPLED_OSC`](TWO_DELAYED_PULSE_COUPLED_OSC/) iterates a
  two-oscillator phase map with delayed pulses.
- [`THREE_DELAYED_PULSE_COUPLED_OSC`](THREE_DELAYED_PULSE_COUPLED_OSC/) extends
  delayed pulse coupling to three oscillators and plots their event times.
- [`TWO_THETA_NEURONS`](TWO_THETA_NEURONS/) classifies a grid of delayed,
  pulse-coupled theta-neuron simulations as synchronized or unsynchronized and
  plots that region in the $(\epsilon,\delta)$ parameter plane.

## What to look for

Vary the initial conditions mentally while following the event sequence: a
locked state repeats the same pulse-arrival phases. In the three-cell plot,
distinguish a repeating collective order from exact simultaneous spikes. For
the theta pair, read the red and blue grid points as synchronized and
unsynchronized outcomes, respectively, and compare them with the plotted
boundary in the $(\epsilon,\delta)$ plane.

## Suggested order

1. Run `TWO_DELAYED_PULSE_COUPLED_OSC`.
2. Run `THREE_DELAYED_PULSE_COUPLED_OSC` and identify its repeating order.
3. Run `TWO_THETA_NEURONS` to connect the phase-map result to continuous
   neuron dynamics.

## Prerequisites and related chapters

Chapter 26 establishes undelayed pulse-coupled phase maps, and Chapter 8
introduces theta neurons. Chapter 28 treats weak coupling, while Chapter 29
focuses on the stability of a synchronous phase relation.

## Running the examples

Run `python main.py` from each immediate example directory. The scripts use
NumPy and Matplotlib; the theta simulation runs a long fine-step trajectory to
classify each sampled coupling-strength and delay pair.
