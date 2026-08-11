# Entrainment by excitatory input pulses

## Overview

Periodic excitatory pulses can make a neuron's spikes adopt the drive's timing.
The examples compare LIF and WB neurons, display one-to-one and n-to-one
responses, and use phase return maps to distinguish locking from irregular
responses.

## Core ideas

Entrainment is phase locking to an external period $T$. In one-to-one
locking, one spike is associated with each input cycle; n-to-one locking has n
input cycles per neuronal cycle or, depending on the event convention, a
repeating multi-cycle relationship. A return map takes the phase after one
pulse to the phase after the next; a stable fixed point predicts a repeatable
phase, while a nonconvergent orbit predicts irregular timing.

## Essential model

If $\alpha_k$ is the phase of the kth relevant spike relative to the pulse
train, the pulse response gives a map

$$
\alpha_{k+1}=F(\alpha_k)\pmod 1.
$$

The detailed WB simulations implement the pulse through a synaptic gate, while
the compact map scripts plot $F$ directly.

## Code examples

All seven examples now live in one notebook, [`chapter23.ipynb`](chapter23.ipynb):
`simulate_lif_entrainment` shows periodic pulses entraining an LIF neuron.
`plot_f_entrainment` plots a phase return map and its intersections with the
identity line. `simulate_f_entrainment_2`/`plot_f_entrainment_2` iterate a
second return-map construction from a chosen initial phase.
`simulate_wb_entrainment_intervals` measures WB spike intervals under
periodic excitation, sweeping the synaptic strength; its inner loop is
JIT-compiled with numba. `simulate_wb_neuron_entrained` gives WB trajectories
with a stable locked response and phases relative to the pulse period.
`simulate_wb_neuron_irregular` demonstrates a WB response whose
pulse-relative timing does not settle into the same pattern.
`simulate_wb_neuron_n_to_one` illustrates an n-to-one WB entrainment pattern.

## What to look for

On a return map, a fixed point lies where $F(\alpha)=\alpha$; repeated
iterations should approach it only when it is stable. In the WB plots, compare
the sequence of pulse-relative spike phases: a locked response repeats, whereas
the irregular trace drifts. For n-to-one locking, count pulse intervals rather
than assuming every pulse evokes a spike.

## Suggested order

1. Run `simulate_lif_entrainment`, `plot_f_entrainment`, and
   `plot_f_entrainment_2`.
2. Run `simulate_wb_neuron_entrained` and `simulate_wb_entrainment_intervals`.
3. Contrast `simulate_wb_neuron_irregular` with `simulate_wb_neuron_n_to_one`.

## Prerequisites and related chapters

Chapter 20 supplies the synaptic-gate mechanism and Chapter 7 supplies the
LIF reset model. Chapters 25--27 make the same phase-map viewpoint explicit
for PRCs, coupled oscillators, and delays.

## Running the examples

Open [`chapter23.ipynb`](chapter23.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom. The
`simulate_wb_entrainment_intervals` cell is JIT-compiled with numba, so
after the first (one-time compile) call it takes well under a minute
instead of the uncompiled sweep's roughly an hour.
