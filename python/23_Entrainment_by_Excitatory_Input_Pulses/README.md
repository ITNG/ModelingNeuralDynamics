# Entrainment by excitatory input pulses

## Overview

Periodic excitatory pulses can make a neuron's spikes adopt the drive's timing.
The examples compare LIF and WB neurons, display one-to-one and n-to-one
responses, and use phase return maps to distinguish locking from irregular
responses.

## Core ideas

Entrainment is phase locking to an external period \(T\). In one-to-one
locking, one spike is associated with each input cycle; n-to-one locking has n
input cycles per neuronal cycle or, depending on the event convention, a
repeating multi-cycle relationship. A return map takes the phase after one
pulse to the phase after the next; a stable fixed point predicts a repeatable
phase, while a nonconvergent orbit predicts irregular timing.

## Essential model

If \(\alpha_k\) is the phase of the kth relevant spike relative to the pulse
train, the pulse response gives a map

\[
\alpha_{k+1}=F(\alpha_k)\pmod 1.
\]

The detailed WB simulations implement the pulse through a synaptic gate, while
the compact map scripts plot \(F\) directly.

## Code examples

- [`LIF_ENTRAINMENT`](LIF_ENTRAINMENT/) shows periodic pulses entraining an LIF
  neuron.
- [`PLOT_F_ENTRAINMENT`](PLOT_F_ENTRAINMENT/) plots a phase return map and its
  intersections with the identity line.
- [`PLOT_F_ENTRAINMENT_2`](PLOT_F_ENTRAINMENT_2/) iterates a second return-map
  construction from a chosen initial phase.
- [`WB_ENTRAINMENT_INTERVALS`](WB_ENTRAINMENT_INTERVALS/) measures WB spike
  intervals under periodic excitation.
- [`WB_NEURON_ENTRAINED`](WB_NEURON_ENTRAINED/) gives a WB trajectory with a
  stable locked response and phases relative to the pulse period.
- [`WB_NEURON_IRREGULAR`](WB_NEURON_IRREGULAR/) demonstrates a WB response
  whose pulse-relative timing does not settle into the same pattern.
- [`WB_NEURON_N_TO_ONE`](WB_NEURON_N_TO_ONE/) illustrates an n-to-one WB
  entrainment pattern.

## What to look for

On a return map, a fixed point lies where \(F(\alpha)=\alpha\); repeated
iterations should approach it only when it is stable. In the WB plots, compare
the sequence of pulse-relative spike phases: a locked response repeats, whereas
the irregular trace drifts. For n-to-one locking, count pulse intervals rather
than assuming every pulse evokes a spike.

## Suggested order

1. Run `LIF_ENTRAINMENT`, `PLOT_F_ENTRAINMENT`, and `PLOT_F_ENTRAINMENT_2`.
2. Run `WB_NEURON_ENTRAINED` and `WB_ENTRAINMENT_INTERVALS`.
3. Contrast `WB_NEURON_IRREGULAR` with `WB_NEURON_N_TO_ONE`.

## Prerequisites and related chapters

Chapter 20 supplies the synaptic-gate mechanism and Chapter 7 supplies the
LIF reset model. Chapters 25--27 make the same phase-map viewpoint explicit
for PRCs, coupled oscillators, and delays.

## Running the examples

Run `python main.py` from the chosen immediate example directory. The WB
simulations use NumPy and Matplotlib; long interval runs use a fine time step
and take longer than the analytical map plots.
