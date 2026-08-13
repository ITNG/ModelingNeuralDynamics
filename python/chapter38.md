# Gamma coherence

## Overview

Gamma coherence is the alignment of population responses across pulses or gamma
episodes. These examples compare two-cell timing references with Poisson-PING
phase sweeps that measure the response of E-cells 1-5, plus a separate raster
with a periodic green waveform as a timing reference.

## Core ideas

Coherence is population timing alignment, not simply high firing rate. A pulse
at a compatible gamma phase can reinforce a common response, while a
mismatched pulse can land during inhibition or at a conflicting phase. Poisson
PING makes this distinction clear because individual spikes vary while a
phase-dependent population response can remain.

## Essential model

For each sampled pulse phase $\varphi$, the Poisson-PING sweep functions
simulate a population, count spikes from E-cells 1-5, and fit a line to the
phase-versus-count response. The two sweep examples use different pulse
periods; the green-reference example instead displays a Poisson-PING raster
with its periodic waveform.

## Code examples

All five examples live in one notebook, [`chapter38.ipynb`](chapter38.ipynb).
A shared model block defines the RTM/WB gating functions (`m_e_inf`,
`h_e_inf`, `tau_h_e`, `n_e_inf`, `tau_n_e`, and their `_i` inhibitory-cell
counterparts) and the synaptic-delay solver `tau_d_q`, reused by every
example below.

- `simulate_two_cell`/`simulate_gamma_coherence_1`/`plot_gamma_coherence_1`
  plot two-cell inhibitory currents with E/I spike times and periodic timing
  markers.
- `simulate_gamma_coherence_2`/`plot_gamma_coherence_2` compare coupled
  inhibition with its mean-inhibition approximation (reusing
  `simulate_two_cell`).
- `simulate_poisson_ping`/`_poisson_ping_loop` (numba-accelerated) build the
  Chapter 38 all-to-all Poisson-PING network, shared by all three
  `POISSON_PING_3_*` examples below.
- `simulate_poisson_ping_mismatched_pulses` sweeps pulse phase for a 29 ms
  period (`P_MISMATCHED_PULSES`) via the shared
  `simulate_poisson_ping_phase_sweep`, plotting E-cell-1-5 spike counts with
  a linear fit.
- `simulate_poisson_ping_plus_green`/`plot_poisson_ping_plus_green` display a
  Poisson-PING E/I raster with a green periodic timing waveform.
- `simulate_poisson_ping_plus_pulses` sweeps pulse phase for a 31 ms period
  (`P_PLUS_PULSES`), reusing the same shared sweep helper.

## What to look for

For the sweep plots, compare E-cell-1-5 spike counts across phase and inspect
the fitted slope rather than reading them as rasters. Use the green waveform
and raster only in its dedicated example to relate a population episode to
the periodic timing reference.

## Suggested order

1. Run `simulate_gamma_coherence_1`/`plot_gamma_coherence_1` and
   `simulate_gamma_coherence_2`/`plot_gamma_coherence_2`.
2. Run `simulate_poisson_ping_plus_green`/`plot_poisson_ping_plus_green` to
   inspect its raster and waveform.
3. Compare the two phase-versus-spike-count sweeps
   (`simulate_poisson_ping_mismatched_pulses` and
   `simulate_poisson_ping_plus_pulses`).

## Prerequisites and related chapters

Chapter 23 introduces pulse timing, Chapter 30 PING, Chapter 32 Poisson PING,
and Chapter 37 thresholding.

## Running the examples

Open [`chapter38.ipynb`](chapter38.ipynb) in Jupyter, or via the Colab badge
at the top of the notebook, and run all cells top to bottom. Poisson-PING
output varies across runs unless a `seed` is fixed, so compare qualitative
timing rather than exact spike identities.
