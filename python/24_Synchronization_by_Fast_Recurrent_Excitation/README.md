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

\[
I_{{\rm syn},i}=g_{\rm syn}\sum_{j\ne i}s_j(v_{\rm syn}-v_i).
\]

The scripts build phase-spread initial conditions from the single-cell limit
cycle and then integrate the full conductance-based network.

## Code examples

- [`RTM_E_TO_E_NETWORK_1`](RTM_E_TO_E_NETWORK_1/) simulates a baseline
  recurrent E-to-E RTM network with release and synaptic-gate variables.
- [`RTM_E_TO_E_NETWORK_2`](RTM_E_TO_E_NETWORK_2/) simulates a longer fast
  E-to-E network run initialized with phases spread across the cycle.
- [`RTM_E_TO_E_HETEROGENEOUS`](RTM_E_TO_E_HETEROGENEOUS/) tests synchrony when
  the network is not identical.
- [`RTM_TWO_CELL_NETWORK`](RTM_TWO_CELL_NETWORK/) makes the phase interaction
  visible in a two-cell reciprocally excitatory network.
- [`RTM_SYNC`](RTM_SYNC/) shows a network started from identical states, giving
  the synchronous reference raster.
- [`RTM_SPLAY`](RTM_SPLAY/) initializes evenly spaced phases to display the
  asynchronous splay configuration.

## What to look for

Compare the aligned raster of `RTM_SYNC` with the evenly staggered spikes of
`RTM_SPLAY`. In the recurrent network, see whether initially different phases
move together over successive cycles. Then inspect the heterogeneous case for
phase dispersion or a shifted collective timing rather than assuming perfect
coincidence is required.

## Suggested order

1. Run `RTM_SYNC` and `RTM_SPLAY` as reference initial conditions.
2. Run `RTM_TWO_CELL_NETWORK`, then `RTM_E_TO_E_NETWORK_1` and
   `RTM_E_TO_E_NETWORK_2`.
3. Finish with `RTM_E_TO_E_HETEROGENEOUS`.

## Prerequisites and related chapters

Chapter 20 develops the fast chemical synapse used here. Chapter 23 shows how
external pulses entrain one neuron; Chapters 25--29 reduce related network
questions to PRCs, phase maps, weak coupling, and stability.

## Running the examples

Run `python main.py` from an immediate example directory. The large network
simulations use NumPy, SciPy or Matplotlib depending on the script, and the
long recurrent run can take noticeably longer than the two-cell example.
