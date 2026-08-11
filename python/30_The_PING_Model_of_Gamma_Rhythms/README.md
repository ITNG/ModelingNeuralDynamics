# The PING model of gamma rhythms

## Overview

PING is a pyramidal--interneuron gamma rhythm: E cells recruit I cells, then recurrent inhibition ends the E-cell episode until recovery starts the next cycle. The examples go from a two-cell loop to sparse and random populations, with rasters and LFP-like output.

## Core ideas

E-to-I excitation is followed by I-to-E inhibition. Drive strength, E/I heterogeneity, and recurrent E/I connectivity decide which cells participate. Population alignment is read from E and I rasters together with mean-voltage (LFP-like) traces.

## Essential model

The E and I cells receive conductance currents such as $I_{\rm IE}=g_{\rm IE}s_I(E_I-v_E)$ and $I_{\rm EI}=g_{\rm EI}s_E(E_E-v_I)$. Random-network weights are normalized by expected in-degree, allowing sparse and dense networks to be compared.

## Code examples

- [`2_CELL_PING`](2_CELL_PING/) integrates the two-cell E--I mechanism and reports the E period.
- [`2_CELL_PING_CONDITION_NUMBERS`](2_CELL_PING_CONDITION_NUMBERS/) perturbs drive, inhibition, and decay to report period sensitivity.
- [`PING_1`](PING_1/) runs a heterogeneous random population and writes E/I spikes and an LFP trace.
- [`PING_2`](PING_2/) runs the corresponding fully connected population.
- [`PING_3`](PING_3/) uses sparse random connectivity.
- [`PING_4`](PING_4/) scales up the sparse population.
- [`PING_5`](PING_5/), [`PING_6`](PING_6/), [`PING_7`](PING_7/), [`PING_8`](PING_8/), and [`PING_9`](PING_9/) plot further drive-dependent population cases.

## What to look for

E spikes should precede I spikes, followed by a silent inhibitory interval. Compare raster bands with LFP peaks; connectivity and drive should change participation and timing, not simply voltage amplitude. Use the condition numbers to spot parameter-sensitive periods.

## Suggested order

1. Run `2_CELL_PING` and `2_CELL_PING_CONDITION_NUMBERS`.
2. Compare `PING_1`, `PING_2`, and `PING_3`, then run `PING_4`.
3. Compare the plotted drive cases `PING_5` through `PING_9`.

## Prerequisites and related chapters

Chapter 20 introduces chemical synapses; Chapters 23--29 introduce phase and synchrony tools. Chapter 31 contrasts all-inhibitory ING, and Chapter 32 studies weak and stochastic PING.

## Running the examples

Run `python main.py` from an immediate example directory. Population runs may write LFP and spike-time outputs used by raster helpers. NumPy, SciPy, and Matplotlib are required.
