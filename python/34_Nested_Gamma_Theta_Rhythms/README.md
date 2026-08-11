# Nested gamma--theta rhythms

## Overview

These examples place gamma PING episodes inside a slower theta cycle. They first develop pre-OLM and OLM-cell h- and A-current dynamics, then combine E, I, and OLM populations under theta-modulated excitation or inhibition.

## Core ideas

Theta modulation selects windows for gamma. h-current promotes slow recovery, while A-current can delay OLM spiking; together they shape OLM inhibition and the timing of E--I activity. Nesting must be read from population rasters relative to the slow modulation.

## Essential model

The E--I--OLM network keeps conductance currents $gs(E_{\rm rev}-v)$ and adds slow OLM gates for h- and A-currents. A theta-periodic E drive or inhibitory conductance creates gamma packets at preferred theta phases.

## Code examples

- [`A_CURRENT`](A_CURRENT/) plots the A-current contribution.
- [`EIO_1`](EIO_1/) simulates an E--I--OLM network and includes a raster helper.
- [`OLM_WITH_H_AND_A_CURRENTS`](OLM_WITH_H_AND_A_CURRENTS/) integrates both slow OLM currents.
- [`OLM_WITH_H_CURRENT`](OLM_WITH_H_CURRENT/) is the h-current-only comparison.
- [`PING_WITH_THETA_DRIVE`](PING_WITH_THETA_DRIVE/) applies theta-modulated excitation.
- [`PING_WITH_THETA_INHIBITION`](PING_WITH_THETA_INHIBITION/) applies theta-modulated inhibition and writes raster/LFP outputs.
- [`PRE_OLM_VOLTAGE_TRACE`](PRE_OLM_VOLTAGE_TRACE/) shows pre-OLM voltage.
- [`PRE_OLM_X_INF_TAU_X`](PRE_OLM_X_INF_TAU_X/) plots pre-OLM gate steady states and time constants.

## What to look for

Compare h-only with h-plus-A OLM dynamics first. Then locate gamma raster packets relative to theta and compare drive-selected with inhibition-selected windows. Use EIO output to separate E, I, and OLM roles.

## Suggested order

1. Run the two `PRE_OLM_*` examples and `A_CURRENT`.
2. Compare the two OLM-current examples.
3. Run `EIO_1`, `PING_WITH_THETA_DRIVE`, and `PING_WITH_THETA_INHIBITION`.

## Prerequisites and related chapters

Chapter 9 introduces slow conductances, Chapter 20 chemical synapses, and Chapters 30--33 PING, ING, weak PING, and beta timing.

## Running the examples

Run `python main.py` from an immediate example directory. Some network cases write LFP and spike-time output for bundled raster helpers. NumPy, SciPy, and Matplotlib are required.
