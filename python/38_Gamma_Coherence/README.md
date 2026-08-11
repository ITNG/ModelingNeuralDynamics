# Gamma coherence

## Overview

Gamma coherence is the alignment of population responses across pulses or gamma
episodes. These examples compare two-cell timing references with Poisson-PING
phase sweeps that measure the response of E-cells 1--5, plus a separate raster
with a periodic green waveform as a timing reference.

## Core ideas

Coherence is population timing alignment, not simply high firing rate. A pulse
at a compatible gamma phase can reinforce a common response, while a
mismatched pulse can land during inhibition or at a conflicting phase. Poisson
PING makes this distinction clear because individual spikes vary while a
phase-dependent population response can remain.

## Essential model

For each sampled pulse phase $\varphi$, the Poisson-PING sweep scripts
simulate a population, count spikes from E-cells 1--5, and fit a line to the
phase-versus-count response. The two sweep directories use different pulse
periods; the green-reference directory instead displays a Poisson-PING raster
with its periodic waveform.

## Code examples

- [`GAMMA_COHERENCE_1`](GAMMA_COHERENCE_1/) plots two-cell inhibitory
  currents with E/I spike times and periodic timing markers.
- [`GAMMA_COHERENCE_2`](GAMMA_COHERENCE_2/) compares coupled inhibition with
  its mean-inhibition approximation.
- [`POISSON_PING_3_MISMATCHED_PULSES`](POISSON_PING_3_MISMATCHED_PULSES/)
  sweeps pulse phase for a 29 ms period and plots E-cell-1--5 spike counts
  with a linear fit.
- [`POISSON_PING_3_PLUS_GREEN`](POISSON_PING_3_PLUS_GREEN/) displays a
  Poisson-PING E/I raster with a green periodic timing waveform.
- [`POISSON_PING_3_PLUS_PULSES`](POISSON_PING_3_PLUS_PULSES/) sweeps pulse
  phase for a 31 ms period and plots E-cell-1--5 spike counts with a linear
  fit.

## What to look for

For the sweep plots, compare E-cell-1--5 spike counts across phase and inspect
the fitted slope rather than reading them as rasters. Use the green waveform
and raster only in its dedicated directory to relate a population episode to
the periodic timing reference. The __pycache__ directory is deliberately
excluded because it is not an example.

## Suggested order

1. Run `GAMMA_COHERENCE_1` and `GAMMA_COHERENCE_2`.
2. Run `POISSON_PING_3_PLUS_GREEN` to inspect its raster and waveform.
3. Compare the two phase-versus-spike-count sweeps.

## Prerequisites and related chapters

Chapter 23 introduces pulse timing, Chapter 30 PING, Chapter 32 Poisson PING, and Chapter 37 thresholding.

## Running the examples

Run `python main.py` in an immediate example directory. NumPy and Matplotlib are required. Poisson-PING output varies across runs, so compare qualitative timing rather than exact spike identities.
