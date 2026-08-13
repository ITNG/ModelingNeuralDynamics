# M-current PING and PINB

## Overview

Chapter 33 covers the book's beta-rhythm material. The examples examine
M-current PING period skipping and PING-with-inhibition/burst (PINB),
including gap junctions and cell-assembly timing.

## Core ideas

Slow M-current recovery can make E cells skip gamma opportunities,
yielding beta-scale population timing. PINB organizes which E assemblies
participate and when; gap junctions can change inhibitory coordination. A
slower rhythm therefore need not mean that every cell fires slowly.

## Essential model

The PING E-I loop is supplemented by M-current-controlled E-cell
availability:

$$
I_M = \hat g_M\, w\,(v_K - v), \qquad
\dot w = \frac{w_\infty(v) - w}{\tau_w(v)}.
$$

Population period and skipped cycles are useful timing measures.
Gap-junction currents are proportional to $v_j-v_i$, providing electrical
coupling without chemical-synapse delay. PINB networks additionally split
the I-cell output into two independently-gated synapse pools -- a fast
one onto other I-cells and a (typically slower or stronger) one onto
E-cells -- both driven by the same I-cell spikes.

## Code examples

All nine examples live in one notebook, [`chapter33.ipynb`](chapter33.ipynb).
The network sims integrate 40-250 coupled cells over tens of thousands of
time steps, so their per-timestep update loops are `@njit`-compiled
(numba), while population initializers and the `tau_d_q` bisection search
stay plain NumPy.

`simulate_m_current_beta_with_gj`/`plot_beta_gj_raster` cover
`M_CURRENT_BETA_WITH_GJ`: a gap-junction-coupled population of M-current E
cells with no synaptic coupling, exposing `g_hat_gap` through an
`interact()` slider.

`simulate_m_current_ping` is the shared M-current PING network integrator
(RTM E cells with M-current, WB I cells, configurable E-E/E-I/I-E/I-I
connectivity). `simulate_m_current_ping_4` through `simulate_m_current_ping_8`
are thin wrappers around it, each exposing one differentiating parameter
through `interact()`: `g_m` (`M_CURRENT_PING_4`), `g_hat_ee`
(`M_CURRENT_PING_5` and `M_CURRENT_PING_6`, the latter with recurrent
excitation confined to the E_P sub-assembly), `i_ext_i_mean`
(`M_CURRENT_PING_7`, with E_S-to-I synapses halved), and the E_P
drive-boost factor (`M_CURRENT_PING_8`). All are plotted with
`plot_m_current_ping_raster`.

`simulate_pinb` is the shared PINB network integrator (plain RTM/WB cells,
no M-current, with independent I-I and I-E synapse pools driven by the
same I-cell spikes). `simulate_pinb_1` (slider: `tau_d_ie`),
`simulate_pinb_2` (slider: `g_hat_ie`), and `simulate_pinb_3` (slider:
`i_ext_e_mean`) are thin wrappers, plotted with `plot_pinb` (raster plus
the E-cell-averaged LFP-like trace).

## What to look for

Count skipped gamma opportunities before calling the population rhythm
beta. In the gap-junction case, inspect cluster timing, not only voltage.
In PINB rasters, identify participating E assemblies and their order
relative to inhibition.

## Suggested order

1. Run `M_CURRENT_PING_4` through `M_CURRENT_PING_8`.
2. Compare `M_CURRENT_BETA_WITH_GJ`.
3. Run `PINB_1`, `PINB_2`, and `PINB_3`.

## Prerequisites and related chapters

Chapter 9 introduces M-current context, Chapter 21 gap junctions, Chapters
30 and 32 PING, and Chapter 34 nested slow/fast rhythms.

## Running the examples

Open [`chapter33.ipynb`](chapter33.ipynb) and run the cells top to bottom.
The network sims (`M_CURRENT_PING_4`-`M_CURRENT_PING_8` and
`PINB_1`-`PINB_3`) integrate 40-250 coupled cells over 200-500ms with a
fixed 0.01ms step, so each can take from several seconds to about a
minute once numba has compiled the step loop.
