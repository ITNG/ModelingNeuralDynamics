# Nested gamma-theta rhythms

## Overview

These examples place gamma PING episodes inside a slower theta cycle. They
first develop pre-O-LM and O-LM cell h- and A-current dynamics, then combine
E, I, and O-LM populations under theta-modulated excitation or inhibition.

## Core ideas

Theta modulation selects windows for gamma. The O-LM cell's h-current
promotes slow post-inhibitory recovery, while its A-current can delay OLM
spiking; together they shape O-LM inhibition and the timing of E-I activity.
Nesting must be read from population rasters relative to the slower
modulation, whether it arrives as a periodic boost to E-cell drive or a
periodic pulse of external inhibition.

## Essential model

The E-I-O-LM network keeps conductance currents $g\,s\,(E_{\rm rev}-v)$ and
adds slow O-LM gates for h- and A-currents. A theta-periodic E drive or
inhibitory conductance creates gamma packets at preferred theta phases.

## Code examples

All eight examples live in one notebook, [`chapter34.ipynb`](chapter34.ipynb).
A shared block of pre-O-LM/O-LM gating functions (`alpha_h`/`alpha_m`/
`alpha_n`/`beta_h`/`beta_m`/`beta_n`, `h_inf`/`m_inf`/`n_inf`,
`tau_h`/`tau_m`/`tau_n`, plus the h-current pair `r_inf`/`tau_r` and the
A-current pair/quad `a_inf`/`b_inf`/`tau_a`/`tau_b`), the RTM/WB E-/I-cell
gates (`m_e_inf`/`h_e_inf`/`tau_h_e`/`n_e_inf`/`tau_n_e`,
`m_i_inf`/`h_i_inf`/`tau_h_i`/`n_i_inf`/`tau_n_i`), the synapse rise-time
solver `tau_peak_function`/`tau_d_q_function`, and the population
splay-state initializers `rtm_init_population`/`wb_init_population`/
`olm_init_population` are defined once and reused throughout the chapter:

- `plot_a_current_gates` plots the O-LM cell's A-current gates and time
  constants (A-current).
- `simulate_eio_network`/`plot_eio_raster` simulate a small E-I-O-LM
  population network (`odeint`-integrated, matching the original script's
  synaptic-term sign convention) and draw its three-population raster plus
  mean($v_E$)/mean($s_E$) traces (E-I-O-LM network, EIO_1).
- `simulate_olm_h_and_a_currents`/`plot_olm_h_and_a_currents` integrate a
  single O-LM cell with both slow currents active (O-LM cell with h- and
  A-currents).
- `simulate_olm_h_current`/`plot_olm_h_current` integrate the h-current
  -only comparison (O-LM cell with h-current).
- `simulate_ping_theta_drive`/`plot_ping_theta_raster` apply a
  theta-modulated multiplicative boost to the E-cell drive (PING with
  theta-modulated drive), with the network stepper numba-accelerated.
- `simulate_ping_theta_inhibition` (sharing `plot_ping_theta_raster`)
  applies a theta-periodic external inhibitory conductance pulse to the E
  cells instead (PING with theta-modulated inhibition).
- `simulate_pre_olm_voltage_trace`/`plot_pre_olm_voltage_trace` show a
  single pre-O-LM cell's voltage trace under constant drive (pre-O-LM
  voltage trace).
- `plot_pre_olm_gates` plots the pre-O-LM cell's steady-state gates and time
  constants (pre-O-LM gating variables).

Interactive sliders let you sweep the O-cell$\to$E-cell coupling strength
`g_hat_oe` in the E-I-O-LM network, the A-/h-current strengths `g_A`/`g_h`
in the two single-cell O-LM examples, the theta modulation depth `alpha` and
external-inhibition strength `g_ex` in the two PING+theta examples, and the
drive `i_ext` in the pre-O-LM voltage trace.

## What to look for

Compare h-only with h-plus-A O-LM dynamics first. Then locate gamma raster
packets relative to theta and compare drive-selected with
inhibition-selected windows. Use the E-I-O-LM network's raster to separate
E, I, and O-LM roles.

## Suggested order

1. Run the pre-O-LM gating/voltage-trace examples and the A-current plot.
2. Compare the two O-LM-current examples (h-only vs. h-plus-A).
3. Run the E-I-O-LM network, then the theta-drive and theta-inhibition PING
   examples.

## Prerequisites and related chapters

Chapter 9 introduces slow conductances, Chapter 20 chemical synapses, and
Chapters 30-33 PING, ING, weak PING, and beta timing.

## Running the examples

Open [`chapter34.ipynb`](chapter34.ipynb) in Jupyter, or via the Colab badge
at the top of the notebook, and run all cells top to bottom. The E-I-O-LM
network (`simulate_eio_network`) integrates a 115-dimensional system with
`scipy.integrate.odeint` and is the slowest cell in the chapter, taking a
few minutes; the theta-modulated PING network cells are numba-accelerated
and much faster.
