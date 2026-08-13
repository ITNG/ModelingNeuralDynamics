# The PING model of gamma rhythms

## Overview

PING is a pyramidal-interneuron gamma rhythm: E cells recruit I cells, then
recurrent inhibition ends the E-cell episode until recovery starts the next
cycle. The examples go from a two-cell loop to sparse and random populations,
with rasters and LFP-like output.

## Core ideas

E-to-I excitation is followed by I-to-E inhibition. Drive strength, E/I
heterogeneity, and recurrent E/I connectivity decide which cells participate.
Population alignment is read from E and I rasters together with mean-voltage
(LFP-like) traces.

## Essential model

The E and I cells receive conductance currents such as
$I_{\rm IE}=g_{\rm IE}s_I(E_I-v_E)$ and $I_{\rm EI}=g_{\rm EI}s_E(E_E-v_I)$.
Random-network weights are normalized by expected in-degree, allowing sparse
and dense networks to be compared.

## Code examples

All eleven examples live in one notebook, [`chapter30.ipynb`](chapter30.ipynb).
A shared block of RTM/WB gating functions, `tau_peak_function`/
`tau_d_q_function` (synapse rise-time solver), `rtm_init_population`/
`wb_init_population` (population splay-state initializers), and
`make_random_connectivity`/`make_fixed_degree_connectivity` (connectivity
builders) is defined once and reused throughout the chapter, along with the
numba-accelerated network stepper `simulate_ping_network` (and its plotting
companions `plot_ping_drive`/`plot_ping_panels`), which is the shared core of
every population example (PING_1-PING_9):

- `simulate_2_cell_ping`/`plot_2_cell_ping` integrate the two-cell E-I
  mechanism and report the E period (2-Cell PING).
- `run_condition_numbers` perturbs drive, inhibition, and decay to report
  period sensitivity (2-Cell PING Condition Numbers).
- `simulate_ping_population` (PING_1-PING_4) runs a random E/I population and
  returns E/I spike times plus an LFP-like trace; called with different
  heterogeneity (`sigma_e`), connectivity (`p_XY`), and size (`num_e`,
  `num_i`) for each of PING_1 (heterogeneous, 50% connectivity), PING_2
  (homogeneous, all-to-all), PING_3 (heterogeneous, sparse), and PING_4
  (PING_3 scaled to 800/200 cells).
- `run_ping5_panels` (PING_5) and `run_connectivity_panels`/`main` (PING_6,
  a plain-NumPy reference version via `simulate_ping_network_plain`) each
  compare dense random, sparse random, and sparse fixed-degree connectivity.
- `simulate_ping_drive` (PING_7, PING_8, PING_9) is a generalized wrapper
  exposing I-cell drive strength, E-E coupling, and I-cell initialization
  mode as parameters.

Interactive sliders let you sweep the I-E coupling strength in the 2-cell
loop, the connectivity density `p` in the PING_1/PING_3 population, and the
mean I-cell drive in the PING_7/PING_8 population.

## What to look for

E spikes should precede I spikes, followed by a silent inhibitory interval.
Compare raster bands with LFP peaks; connectivity and drive should change
participation and timing, not simply voltage amplitude. Use the condition
numbers to spot parameter-sensitive periods.

## Suggested order

1. Run the 2-Cell PING and 2-Cell PING Condition Numbers sections.
2. Compare PING_1, PING_2, and PING_3, then run PING_4.
3. Compare the plotted drive cases PING_5 through PING_9.

## Prerequisites and related chapters

Chapter 20 introduces chemical synapses; Chapters 23-29 introduce phase and
synchrony tools. Chapter 31 contrasts all-inhibitory ING, and Chapter 32
studies weak and stochastic PING.

## Running the examples

Open [`chapter30.ipynb`](chapter30.ipynb) in Jupyter, or via the Colab badge
at the top of the notebook, and run all cells top to bottom. These are
spiking-network simulations and some cells (particularly PING_4 and PING_6's
default 2-second run) can take a while to finish -- that is expected.
