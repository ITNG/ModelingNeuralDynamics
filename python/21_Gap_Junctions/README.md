# Gap junctions

## Overview

Gap junctions are electrical synapses: each cell receives a current determined
by the voltage difference from its neighbours. The examples contrast the
resulting synchronization of spiking cells with subthreshold coupling. They
also contrast LIF event handling with a continuous WB voltage trace.

## Core ideas

For a pair, a gap-junction current is proportional to \(v_j-v_i\), so it
reduces voltage differences rather than imposing a fixed chemical reversal
potential. This promotes synchrony when cells spike, but subthreshold behavior
can still be distinctive. In LIF models, a spike reset makes the threshold and
reset separation part of the coupling mechanism.

## Essential model

With electrical conductance \(g_{\rm gap}\), cell \(i\) receives

\[
I_{{\rm gap},i}=g_{\rm gap}\sum_{j\ne i}(v_j-v_i).
\]

The continuous current is supplemented by the model-specific event rule: an
LIF voltage is reset after crossing threshold, whereas the WB network evolves
its conductance-based voltage continuously.

## Code examples

- [`LIF_NETWORK_WITH_GJ`](LIF_NETWORK_WITH_GJ/) compares two LIF voltage
  traces under diffusive gap-junction coupling, with and without the script's
  additional spike-triggered voltage kick to the other cell; vertical marks
  identify spikes.
- [`RESET_THRESHOLD`](RESET_THRESHOLD/) integrates a single WB
  conductance-based voltage trace and marks two reference voltage levels.
- [`WB_NETWORK_WITH_GJ`](WB_NETWORK_WITH_GJ/) integrates two WB neurons with a
  gap junction to display their voltage alignment.
- [`WB_NETWORK_WITH_GJ_SUBTHRESHOLD`](WB_NETWORK_WITH_GJ_SUBTHRESHOLD/) focuses
  on voltage-difference coupling before the WB cells spike.

## What to look for

Inspect whether an initial voltage difference shrinks in the subthreshold WB
trace and whether spike times become aligned in the network plots. In the LIF
comparison, separate diffusive electrical equalization from the additional
spike-triggered kick. Relate the discontinuous reset to the voltage gap
immediately after a threshold crossing in `LIF_NETWORK_WITH_GJ`; contrast this
with the smooth WB trajectories, including `RESET_THRESHOLD`.

## Suggested order

1. Run `LIF_NETWORK_WITH_GJ`, then `RESET_THRESHOLD` as a continuous WB
   voltage reference.
2. Run `WB_NETWORK_WITH_GJ_SUBTHRESHOLD` to isolate electrical equalization.
3. Run `WB_NETWORK_WITH_GJ` and compare the spiking case.

## Prerequisites and related chapters

Chapter 7 defines LIF threshold-and-reset dynamics, and the WB examples use
the conductance-based interneuron model from earlier chapters. Chapter 20
covers chemical synapses; Chapters 24 and 29 return to synchrony and its
stability.

## Running the examples

Run `python main.py` from each immediate example directory. NumPy, SciPy, and
Matplotlib are required, and each script writes its figure in that directory.
