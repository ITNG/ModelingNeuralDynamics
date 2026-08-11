# Periodic inhibition

## Overview

Periodic inhibitory forcing creates recurring response windows rather than merely lowering firing. The examples show oscillatory responses, three forcing constructions, and f--I curves under periodic inhibition for simple and RTM settings.

## Core ideas

Spikes can be locked, suppressed, or shifted according to their phase in the inhibitory cycle. Inhibited f--I curves count spikes that survive these windows, so changes in average rate can reflect cycle skipping rather than a uniformly slower intrinsic response.

## Essential model

The forcing is $I_{\rm inh}(t)=g_{\rm inh}s(t)(E_I-v)$, with a periodic gate $s(t)$. Sweeping applied current and counting spikes per observation time produces the inhibited f--I relation.

## Code examples

- [`OSCILLATIONS`](OSCILLATIONS/) plots a reference oscillatory response.
- [`PERIODIC_INHIBITION`](PERIODIC_INHIBITION/), [`PERIODIC_INHIBITION_2`](PERIODIC_INHIBITION_2/), and [`PERIODIC_INHIBITION_3`](PERIODIC_INHIBITION_3/) simulate periodic inhibitory forcing variants.
- [`PERIODIC_INHIBITION_F_I_CURVE`](PERIODIC_INHIBITION_F_I_CURVE/) and [`PERIODIC_INHIBITION_F_I_CURVE_2`](PERIODIC_INHIBITION_F_I_CURVE_2/) compute inhibited f--I curves.
- [`RTM_F_I_CURVE_WITH_INHIBITION`](RTM_F_I_CURVE_WITH_INHIBITION/) and [`RTM_F_I_CURVE_WITH_INHIBITION_2`](RTM_F_I_CURVE_WITH_INHIBITION_2/) give RTM comparisons.

## What to look for

Follow spikes relative to inhibitory cycles, especially at a response-window boundary. Compare f--I slopes, offsets, and plateaus: lower rate can arise from skipped cycles.

## Suggested order

1. Run `OSCILLATIONS` and the three periodic-inhibition traces.
2. Compare the two periodic-inhibition f--I curves.
3. Contrast the two RTM f--I calculations.

## Prerequisites and related chapters

Chapter 17 introduces steady f--I curves, Chapter 23 periodic forcing, Chapter 31 ING, and Chapter 36 pulsed excitation.

## Running the examples

Run `python main.py` in an immediate example directory. NumPy and Matplotlib are required; f--I sweeps take longer than one trace.
