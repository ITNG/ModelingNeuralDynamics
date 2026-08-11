# Synchronization by fast recurrent excitation

## Overview

Fast recurrent excitation can bring initially dispersed RTM neurons into a
common firing rhythm. The examples progress from two-cell and recurrent-network
cases to synchronous, splay-state, and heterogeneous networks.

## Core ideas

An excitatory pulse advances cells that are susceptible to it, so repeated
interactions can contract phase differences. Synchrony means the cells fire
together; a splay state distributes their phases around the cycle. Heterogeneous
drive or intrinsic properties work against a common rhythm, testing how robust
the coupling mechanism is.

## Essential model

Each RTM cell receives summed recurrent excitation through synaptic gates:

$$
I_{{\rm syn},i}=g_{\rm syn}\sum_{j\ne i}s_j(v_{\rm syn}-v_i).
$$

The scripts build phase-spread initial conditions from the single-cell limit
cycle and then integrate the full conductance-based network.

## Code examples

All six examples now live in one notebook, [`chapter24.ipynb`](chapter24.ipynb):
`simulate_rtm_e_to_e_network` simulates a baseline recurrent E-to-E RTM
network with release and synaptic-gate variables; `simulate_rtm_e_to_e_network_2`
calls the same function for a much longer run initialized with phases spread
across the cycle, watching the splay state contract toward synchrony -- its
inner loop is JIT-compiled with numba. `simulate_rtm_e_to_e_heterogeneous`
tests synchrony when the network is not identical (random per-neuron drive
and per-pair coupling). `simulate_rtm_two_cell_network` makes the phase
interaction visible in a two-cell reciprocally excitatory network.
`simulate_rtm_sync` shows a network started from identical states, giving
the synchronous reference raster. `simulate_rtm_splay` initializes evenly
spaced phases to display the asynchronous splay configuration. All the
network examples share one `rtm_init` helper that finds the single-cell
RTM limit cycle and interpolates phase-spread initial conditions from it.

## What to look for

Compare the aligned raster of `simulate_rtm_sync` with the evenly staggered
spikes of `simulate_rtm_splay`. In the recurrent network, see whether
initially different phases move together over successive cycles. Then
inspect the heterogeneous case for phase dispersion or a shifted collective
timing rather than assuming perfect coincidence is required.

## Suggested order

1. Run `simulate_rtm_sync` and `simulate_rtm_splay` as reference initial
   conditions.
2. Run `simulate_rtm_two_cell_network`, then `simulate_rtm_e_to_e_network`
   and `simulate_rtm_e_to_e_network_2`.
3. Finish with `simulate_rtm_e_to_e_heterogeneous`.

## Prerequisites and related chapters

Chapter 20 develops the fast chemical synapse used here. Chapter 23 shows how
external pulses entrain one neuron; Chapters 25--29 reduce related network
questions to PRCs, phase maps, weak coupling, and stability.

## Running the examples

Open [`chapter24.ipynb`](chapter24.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom. The
`simulate_rtm_e_to_e_network_2` cell is JIT-compiled with numba, so after
the first (one-time compile) call the long recurrent run takes seconds
instead of minutes.
