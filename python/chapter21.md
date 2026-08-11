# Gap junctions

## Overview

Gap junctions are electrical synapses: each cell receives a current determined
by the voltage difference from its neighbours. The examples contrast the
resulting synchronization of spiking cells with subthreshold coupling. They
also contrast LIF event handling with a continuous WB voltage trace.

## Core ideas

For a pair, a gap-junction current is proportional to $v_j-v_i$, so it
reduces voltage differences rather than imposing a fixed chemical reversal
potential. This promotes synchrony when cells spike, but subthreshold behavior
can still be distinctive. In LIF models, a spike reset makes the threshold and
reset separation part of the coupling mechanism.

## Essential model

With electrical conductance $g_{\rm gap}$, cell $i$ receives

$$
I_{{\rm gap},i}=g_{\rm gap}\sum_{j\ne i}(v_j-v_i).
$$

The continuous current is supplemented by the model-specific event rule: an
LIF voltage is reset after crossing threshold, whereas the WB network evolves
its conductance-based voltage continuously.

## Code examples

All four examples now live in one notebook, [`chapter21.ipynb`](chapter21.ipynb):
`simulate_lif_network_with_gj` compares two LIF voltage traces under
diffusive gap-junction coupling, with and without an additional
spike-triggered voltage kick to the other cell (`epsilon`); vertical marks
identify spikes. `simulate_reset_threshold` integrates a single WB
conductance-based voltage trace and marks two reference voltage levels.
`simulate_wb_network_with_gj` integrates two WB neurons with a gap junction
to display their voltage alignment. `simulate_wb_network_with_gj_subthreshold`
focuses on voltage-difference coupling before the WB cells spike, taking a
`gap_gate(v1)` callable so the same stepper can run the always-on and
subthreshold-only cases.

## What to look for

Inspect whether an initial voltage difference shrinks in the subthreshold WB
trace and whether spike times become aligned in the network plots. In the LIF
comparison, separate diffusive electrical equalization from the additional
spike-triggered kick. Relate the discontinuous reset to the voltage gap
immediately after a threshold crossing in `simulate_lif_network_with_gj`;
contrast this with the smooth WB trajectories, including
`simulate_reset_threshold`.

## Suggested order

1. Run `simulate_lif_network_with_gj`, then `simulate_reset_threshold` as a
   continuous WB voltage reference.
2. Run `simulate_wb_network_with_gj_subthreshold` to isolate electrical
   equalization.
3. Run `simulate_wb_network_with_gj` and compare the spiking case.

## Prerequisites and related chapters

Chapter 7 defines LIF threshold-and-reset dynamics, and the WB examples use
the conductance-based interneuron model from earlier chapters. Chapter 20
covers chemical synapses; Chapters 24 and 29 return to synchrony and its
stability.

## Running the examples

Open [`chapter21.ipynb`](chapter21.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom.
