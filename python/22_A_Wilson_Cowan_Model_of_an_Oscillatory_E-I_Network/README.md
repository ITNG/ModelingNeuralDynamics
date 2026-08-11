# A Wilson--Cowan model of an oscillatory E--I network

## Overview

Wilson--Cowan equations replace individual membrane voltages with firing-rate
variables for excitatory and inhibitory populations. These examples show how
recurrent excitation and delayed inhibitory feedback make a population rhythm,
and how that rhythm appears in rate traces, a phase plane, and a rastergram.

## Core ideas

The variables $E$ and $I$ are population activities filtered on their own
time scales. Nonlinear response functions turn recurrent input into rates.
Excitation raises the E population, inhibition subsequently suppresses it, and
the cycle can repeat. Lowering recurrent excitation changes whether this loop
can sustain a rhythm.

## Essential model

A representative system is

$$
\tau_E\dot E=-E+f(w_{EE}E-w_{IE}I+P_E),\qquad
\tau_I\dot I=-I+g(w_{EI}E-w_{II}I+P_I).
$$

Nullclines are the points where one derivative is zero; their geometry helps
explain the direction of a trajectory and any oscillatory feedback loop.

## Code examples

- [`WILSON_COWAN_E_AND_I`](WILSON_COWAN_E_AND_I/) plots the excitatory and
  inhibitory rate traces over several oscillation cycles.
- [`WILSON_COWAN_LOWERING_W_EE`](WILSON_COWAN_LOWERING_W_EE/) lowers recurrent
  excitation and shows the altered long-time activity.
- [`WILSON_COWAN_PHASE_PLANE`](WILSON_COWAN_PHASE_PLANE/) draws the E--I
  nullclines and the trajectory in the rate phase plane.
- [`WILSON_COWAN_RASTERGRAM`](WILSON_COWAN_RASTERGRAM/) samples the rates into
  a population-style raster representation.

## What to look for

In `WILSON_COWAN_E_AND_I`, compare the timing of E and I peaks rather than
only their heights. In the phase plane, follow the trajectory relative to both
nullclines. Then compare the lowered-$w_{EE}$ trace and the raster: a rate
oscillation is summarized differently by continuous activity and by events.

## Suggested order

1. Run `WILSON_COWAN_E_AND_I`.
2. Use `WILSON_COWAN_PHASE_PLANE` to interpret the rate cycle.
3. Compare `WILSON_COWAN_LOWERING_W_EE` and `WILSON_COWAN_RASTERGRAM`.

## Prerequisites and related chapters

The earlier phase-plane and bifurcation chapters provide the geometric tools
for rate nullclines. Chapter 24 studies synchronization in spiking excitatory
networks, whereas this chapter describes the population-level E--I mechanism.

## Running the examples

Run `python main.py` in an immediate example directory. They use NumPy, SciPy,
and Matplotlib and save figure files alongside the scripts.
