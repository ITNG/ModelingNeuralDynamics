# Weak PING rhythms

## Overview

The repository's historical title, “M Current PING and Poisson PING,” covers the book's weak-PING material. It compares stochastic Poisson drive with deterministic M-current weak PING through phase functions, closeups, and clustering.

## Core ideas

Poisson PING can make irregular individual events align into gamma-population episodes. M-current PING uses slow potassium-mediated recovery to create deterministic weak PING. Phase functions and cluster plots explain recruitment timing rather than merely displaying a period.

## Essential model

The E cell includes $I_M=g_Mw(v-E_K)$ with $\dot w=(w_\infty(v)-w)/\tau_w(v)$. Poisson cases replace fixed drive with random events, so alignment must be judged from population timing rather than an identical raster on every run.

## Code examples

- [`M_CURRENT_PING_1`](M_CURRENT_PING_1/) is the baseline M-current weak-PING network.
- [`M_CURRENT_PING_1_CLOSEUP`](M_CURRENT_PING_1_CLOSEUP/) zooms its within-cycle timing.
- [`M_CURRENT_PING_1_FROM_REST`](M_CURRENT_PING_1_FROM_REST/) starts it from rest.
- [`M_CURRENT_PING_2_CLOSEUP`](M_CURRENT_PING_2_CLOSEUP/) and [`M_CURRENT_PING_3_CLOSEUP`](M_CURRENT_PING_3_CLOSEUP/) compare later M-current cases.
- [`PING_CLUSTERS`](PING_CLUSTERS/) shows temporal PING clusters.
- [`PLOT_PHI`](PLOT_PHI/), [`PLOT_PSI`](PLOT_PSI/), and [`PLOT_PSI_PHI`](PLOT_PSI_PHI/) plot the phase functions used to interpret timing.
- [`POISSON_PING_1`](POISSON_PING_1/), [`POISSON_PING_2`](POISSON_PING_2/), and [`POISSON_PING_3`](POISSON_PING_3/) are stochastic population cases.
- [`POISSON_PING_3_VOLTAGE_TRACE`](POISSON_PING_3_VOLTAGE_TRACE/) zooms the third Poisson case.

## What to look for

Contrast irregular input events with aligned population episodes. In M-current closeups, follow slow recovery before assigning an E--I cycle. Read the phase plots with `PING_CLUSTERS`: phase structure explains clustered timing.

## Suggested order

1. Run the baseline, from-rest, and first closeup M-current cases.
2. Inspect the other M-current closeups, phase plots, and clusters.
3. Compare all three Poisson cases, then their voltage closeup.

## Prerequisites and related chapters

Chapter 9 provides slow-current context, Chapter 30 PING, Chapter 33 beta-rhythm material, and Chapter 38 Poisson-PING coherence.

## Running the examples

Run `python main.py` from an immediate example directory. NumPy, SciPy, and Matplotlib are required. Poisson output is stochastic, so compare qualitative timing across runs.
