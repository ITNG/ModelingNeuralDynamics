# Gamma coherence

## Overview

Gamma coherence is the alignment of population responses across pulses or gamma episodes. These examples compare matched and mismatched pulse timing and perturb a Poisson-PING population to show stronger, shifted, or disrupted alignment.

## Core ideas

Coherence is population timing alignment, not simply high firing rate. A pulse at a compatible gamma phase can reinforce a common response, while a mismatched pulse can land during inhibition or at a conflicting phase. Poisson PING makes this distinction clear because individual spikes vary while collective phase preference can remain.

## Essential model

For phases \(\phi_k\) relative to gamma, a population-alignment summary is \(R=|N^{-1}\sum_k e^{i\phi_k}|\). The examples express the same idea through rasters and voltage responses to differently timed pulses.

## Code examples

- [`GAMMA_COHERENCE_1`](GAMMA_COHERENCE_1/) shows a first gamma-coherence pulse response.
- [`GAMMA_COHERENCE_2`](GAMMA_COHERENCE_2/) compares a second timing condition.
- [`POISSON_PING_3_MISMATCHED_PULSES`](POISSON_PING_3_MISMATCHED_PULSES/) perturbs Poisson PING with mismatched pulses.
- [`POISSON_PING_3_PLUS_GREEN`](POISSON_PING_3_PLUS_GREEN/) adds a marked green timing reference.
- [`POISSON_PING_3_PLUS_PULSES`](POISSON_PING_3_PLUS_PULSES/) overlays applied pulses on the Poisson-PING output.

## What to look for

Compare pulse time with the E--I rhythm, then look for tighter or more dispersed population response. The marked and pulse-overlay cases distinguish true alignment from extra spikes. The __pycache__ directory is deliberately excluded because it is not an example.

## Suggested order

1. Run `GAMMA_COHERENCE_1` and `GAMMA_COHERENCE_2`.
2. Run the two `POISSON_PING_3_PLUS_*` examples.
3. Contrast them with `POISSON_PING_3_MISMATCHED_PULSES`.

## Prerequisites and related chapters

Chapter 23 introduces pulse timing, Chapter 30 PING, Chapter 32 Poisson PING, and Chapter 37 thresholding.

## Running the examples

Run `python main.py` in an immediate example directory. NumPy and Matplotlib are required. Poisson-PING output varies across runs, so compare qualitative timing rather than exact spike identities.
