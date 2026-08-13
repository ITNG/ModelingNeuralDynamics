# ING rhythms

## Overview

Interneuron gamma (ING) is an inhibitory network rhythm in which recovery
from shared inhibition determines the next spike window. These examples
move from a one-cell timing reference and abstract inhibitory pulse-coupling
maps, through full spiking-network simulations of an inhibitory population
(with optional electrical coupling), to ING entraining a population of
excitatory cells.

## Core ideas

ING synchrony is organized by inhibitory recovery rather than recurrent
excitation. Inhibitory phase-reset maps predict stable timing; gap
junctions can alter clustering. The population outputs show whether an
inhibitory rhythm creates repeatable windows for excitatory spikes.

## Essential model

Recurrent inhibition has the form $I_{{\rm II},i}=g_{\rm II}\sum_j s_j(E_I-v_i)$,
optionally supplemented by gap-junction voltage-difference currents. A
pulse map $\phi\mapsto F(\phi)$ reduces inhibitory timing to fixed points
and their slopes.

## Code examples

All sixteen examples live in one notebook, [`chapter31.ipynb`](chapter31.ipynb).
The `ING_1`-`ING_10` and entrainment population sims share `@njit`-compiled
(numba) per-timestep update loops, and most sub-examples expose one natural
parameter (synaptic/gap-junction strength, heterogeneity, connectivity,
drive, ...) through an `interact()` slider so you can sweep it and re-plot
without re-running the whole notebook.
`simulate_1_cell_ing`/`plot_1_cell_ing` give the single self-inhibited-cell
timing reference. `compute_condition_numbers` measures the one-cell
period's sensitivity to `i_ext`, `g_ii`, and `tau_d`.
`simulate_abstract_pulse_coupling_inh`/`plot_abstract_pulse_coupling_inh`
and `simulate_abstract_pulse_coupling_inh_2`/`plot_abstract_pulse_coupling_inh_2`
plot the two inhibitory pulse-coupling maps. `simulate_ing_population` is
the shared WB inhibitory-population network integrator (random or
fixed-indegree synaptic wiring, optional sparse gap junctions);
`simulate_ing_1` through `simulate_ing_10` are thin wrappers around it with
different heterogeneity/connectivity/gap-junction parameters, plotted with
`plot_ing_raster`, `plot_ing_raster_scalebar` (`ING_7`), or
`plot_ing_raster_zoom` (`ING_8`-`ING_10`). `build_entrainment_network` and
`simulate_entrainment` are the shared RTM-E/WB-I network for the
entrainment examples: `simulate_ing_entraining_e_cells` runs the reference
configuration, and `run_drive_panels`/`main` sweep three E-cell drive
levels for the second entrainment figure.

## What to look for

Cells should fire after inhibitory conductance decays enough to open a
recovery window. Compare abstract-map fixed points with persistent raster
timing. In the entrainment examples, identify E spikes relative to the
inhibitory window rather than expecting every E cell to fire every cycle.

## Suggested order

1. Run the one-cell and abstract-map examples (`simulate_1_cell_ing`,
   `compute_condition_numbers`, the two pulse-coupling maps).
2. Compare `simulate_ing_1` through `simulate_ing_10`.
3. Finish with `simulate_ing_entraining_e_cells` and `run_drive_panels`/`main`.

## Prerequisites and related chapters

Chapter 21 introduces gap junctions, Chapters 25-29 phase maps, Chapter 30
PING, and Chapter 35 periodic inhibition.

## Running the examples

Open [`chapter31.ipynb`](chapter31.ipynb) and run the cells top to bottom.
The population network sims (`ING_1`-`ING_10` and the two entrainment
examples) integrate 100-500 coupled cells over 500ms with a fixed 0.01ms
step, so each can take from several seconds to about a minute.
