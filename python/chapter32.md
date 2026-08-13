# Weak PING rhythms

## Overview

This chapter compares two ways a PING population produces a *weak*,
loosely-periodic rhythm instead of a tightly-locked one: an M-current-
mediated recovery in the E cells ("M-current PING"), and a Poisson stream
of independent excitatory events driving each E cell ("Poisson PING"). It
also derives the phase maps used to interpret temporal sub-clustering
within a PING volley.

## Core ideas

Poisson PING lets irregular, uncorrelated individual events align into
gamma-population episodes through recurrent inhibition. M-current PING
instead uses a slow, non-inactivating potassium current so that an E
cell's own recent spiking suppresses its next spike, independent of the
inhibitory volley. Phase maps `psi`/`phi` explain recruitment timing within
a cycle rather than merely displaying a period.

## Essential model

The E cell adds an M-current $I_M=g_M w(v-E_K)$ with
$\dot w=(w_\infty(v)-w)/\tau_w(v)$ to the usual RTM E cell; the I cell is
the usual WB cell. Poisson examples replace part of the fixed drive with
independent random synaptic events per E cell, so alignment must be judged
from population timing rather than an identical raster on every run.

## Code examples

All thirteen examples live in one notebook, [`chapter32.ipynb`](chapter32.ipynb).
`simulate_m_current_ping_plain` is the shared plain-NumPy M-current network
stepper (`M_CURRENT_PING_1`, and `PING_CLUSTERS` via `g_m=0`);
`simulate_m_current_ping_numba` is its numba-accelerated sibling, built on
the compiled `_m_current_run_loop`/`_m_current_settle_loop`, used by
`simulate_m_current_ping_closeup` (`M_CURRENT_PING_1_CLOSEUP`,
`_2_CLOSEUP`, `_3_CLOSEUP` -- identical parameters in the legacy scripts, so
one function serves all three) and `simulate_m_current_ping_1_from_rest`
(`M_CURRENT_PING_1_FROM_REST`, which settles at rest for 200 ms before the
drive turns on). `simulate_ping_clusters` reuses the plain stepper with
`g_m=0` for `PING_CLUSTERS`. `simulate_plot_phi`, `simulate_plot_psi`, and
`simulate_plot_psi_phi` compute the `psi`/`phi` inhibitory phase-reset
maps. `simulate_poisson_ping` is the shared Poisson-PING stepper (no
M-current, an independent Poisson synaptic stream per E cell);
`simulate_poisson_ping_1`, `_2`, `_3` are thin parameter wrappers around it,
and `simulate_poisson_ping_3` doubles as `POISSON_PING_3_VOLTAGE_TRACE`
since it already returns the tracked E cell's voltage trace. Several
sub-examples expose a natural scalar parameter (`g_m`, `f_stoch`, `g_I`)
through an `interact()` slider.

## What to look for

Contrast irregular Poisson input events with aligned population episodes.
In the M-current closeups, follow the slow `w` recovery before assigning an
E-I cycle; `M_CURRENT_PING_1_FROM_REST`'s raster should stay empty for the
first ~200 ms. Read the phase plots alongside `PING_CLUSTERS`: phase
structure explains clustered timing.

## Suggested order

1. Run `M_CURRENT_PING_1`, `_FROM_REST`, and the first closeup.
2. Inspect the remaining M-current closeups, the phase plots, and
   `PING_CLUSTERS`.
3. Compare all three Poisson-PING cases, then the voltage-trace closeup.

## Prerequisites and related chapters

Chapter 9 provides slow-current context, Chapter 30 PING, Chapter 33
beta-rhythm material, and Chapter 38 Poisson-PING coherence.

## Running the examples

Open [`chapter32.ipynb`](chapter32.ipynb) and run the cells top to bottom.
NumPy, SciPy, Matplotlib, numba, and ipywidgets are required. Poisson
output is stochastic, so compare qualitative timing across runs rather than
exact spike times.
