# Thresholding in PING

## Overview

PING can have sharp boundaries between suppression and participation. This
chapter's notebook isolates a non-reset threshold mechanism, runs a
baseline PING network to locate the boundary, magnifies the transition
where a small timing or drive change flips a cell's participation, and
computes the boundary directly for a single cell driven by a periodic
inhibitory conductance.

## Core ideas

Thresholding is a network property: excitation must arrive during the
interval left open by recurrent inhibition. A non-reset reference cell
separates continuous voltage crossing from an imposed reset. Near a
boundary, a cell can miss a cycle or get recruited into the next one from
what looks like a negligible change in drive or timing.

## Essential model

The E-to-I and I-to-E PING conductance loop evolves continuously, while a
voltage threshold selects spike events. In the non-reset comparison,
crossing does not reset the voltage. In the direct threshold construction,
a single RTM cell is driven by a fixed periodic inhibitory conductance
trace $g(t)=e^{\cos^4(\pi t/25)}-1$ (scaled to a target mean $\bar g$); the
onset drive $I_L$ (first spike) and the drive $I_R$ at which the cell locks
onto the full 39 Hz rhythm bracket a thresholding window $w=I_R-I_L$ that
narrows as $\bar g$ grows.

## Code examples

All four examples live in one notebook, [`chapter37.ipynb`](chapter37.ipynb).

- `no_reset_time_constants`/`simulate_no_reset`/`plot_no_reset` (**NO_RESET**)
  build the piecewise-exponential, never-reset voltage trace and its two
  effective time constants (`tau_m`, `tau_m_hat`).
- `simulate_ping_thr` is the shared, `@njit`-compiled RTM-E/WB-I PING
  network stepper (200 E cells with a linearly ramping drive, 50 I cells,
  all-to-all connectivity) reused by both PING examples below.
  `simulate_ping_thr_1`/`plot_ping_thr_raster` (**PING_THR_1**) run it for
  200 ms and raster the whole population, with the E-cell index doubling
  as a drive axis. `simulate_ping_thr_1_zoom`/`plot_ping_thr_zoom`
  (**PING_THR_1_ZOOM**) run it for 500 ms and zoom into E-cells 72-78 to
  inspect the suppression/participation boundary cycle by cycle.
- `firing_rate`, `bisect_threshold`, and `threshold_width_sweep`
  (**THRESHOLDING**) compute the onset drive `I_L`, the full-rate-locking
  drive `I_R`, and the resulting window `w = I_R - I_L` for five values of
  `g_bar`; the per-drive RTM simulation is `@njit`-compiled since dozens of
  bisection evaluations are needed per `g_bar`.

`interact()` sliders are provided for `g_bar` (NO_RESET) and `g_hat_ie`
(PING_THR_1), the two single-parameter knobs each example is built around.

## What to look for

Use `NO_RESET` to distinguish a threshold event from reset dynamics. In
`PING_THR_1_ZOOM`, inspect whether an E cell crosses when inhibition
relaxes enough to recruit it, or remains suppressed cycle after cycle. A
boundary represents changed cycle participation, not merely a small
voltage difference. In `THRESHOLDING`, note how `w` shrinks by roughly a
factor of 2.5-3 for each 0.05 increment of `g_bar`.

## Suggested order

1. Run `NO_RESET` and `THRESHOLDING`.
2. Run `PING_THR_1`.
3. Use `PING_THR_1_ZOOM` to inspect its boundary.

## Prerequisites and related chapters

Chapter 7 discusses reset dynamics, Chapter 30 PING, Chapter 35 inhibitory
windows, and Chapter 38 gamma coherence.

## Running the examples

Open [`chapter37.ipynb`](chapter37.ipynb) and run the cells top to bottom.
The PING network cells (`PING_THR_1`, `PING_THR_1_ZOOM`) integrate 250
cells over 200-500 ms at a fixed 0.01 ms step, and `THRESHOLDING` runs
dozens of 5000 ms bisection evaluations; each of these can take from
several seconds to about a minute.
