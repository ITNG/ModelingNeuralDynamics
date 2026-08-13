# Phase locking of two oscillators

## Overview

Two pulse-coupled oscillators can be studied with an event-to-event phase map
instead of their full voltage trajectories. The abstract examples build and
iterate maps with symmetric or asymmetric PRCs; the RTM example compares the
map ingredients with a conductance-based neuron.

## Core ideas

A locking relationship is a fixed point of the phase map. Its slope determines
local stability: nearby phase differences return only when the map contracts
them. Symmetric PRCs can support symmetric phase relationships, whereas
asymmetric responses shift or remove them. Directly iterating a pair of pulse
oscillators displays the same stable or unstable phase arrangements.

## Essential model

If $g$ is a phase-resetting curve, the examples first form the single-pulse
reset map $f(\phi)=\phi+g(\phi)$. They then account for the other
oscillator's complementary phase with $F(\phi)=f(1-\phi)$, and form the
two-event phase-difference map $G(\phi)=F(F(\phi))$. Locking satisfies

$$
G(\phi_*)=\phi_*.
$$

For a one-dimensional map, the fixed point is locally stable when
$|G'(\phi_*)|<1$. The examples plot the two-event map and identity line so
those intersections can be read directly.

## Code examples

All nine examples now live in one notebook, [`chapter26.ipynb`](chapter26.ipynb).
`abstract_phase_grid` gives the shared $\varphi\in[0,1]$ grid, and
`pulse_map_f`/`pulse_map_bigF`/`pulse_map_bigG` build $f=\varphi+g(\varphi)$,
$F=f(1-\cdot)$, and $G=F\circ F$ from any PRC function `g`; `check_f_monotonic`
prints a warning if $f$ turns out not to be strictly increasing.
`plot_pulse_coupling_full` (a 2x2 figure of $g,f,F,G$) and
`plot_pulse_coupling_g_and_bigG` (just $g$ and $G$) are the two plotting
helpers reused across the five abstract PRCs: `g_pulse_1`
($\varphi^2(1-\varphi)$, asymmetric positive), `g_pulse_2`
($\epsilon\varphi(1-\varphi)^3$), `g_pulse_3` (antisymmetric about
$\varphi=1/2$), `g_pulse_4` (a shifted-locking asymmetric response), and
`g_pulse_5` (the negative counterpart of `g_pulse_1`). `simulate_f_tilde`/
`plot_f_tilde` plot the transformed map $\tilde f$ used to locate phase-map
fixed points geometrically. `simulate_two_pulse_coupled_osc` iterates two
abstract pulse-coupled oscillators event by event from given starting phases
(using PRCs `g_two_pulse_1` and `g_two_pulse_2`), and
`plot_two_pulse_coupled_osc` draws the resulting spike raster. `rtm_init`
finds the RTM limit cycle for splay initialization; `simulate_rtm_g` computes
the RTM phase-resetting curve from a single synaptic pulse (the same recipe
as Chapter 25's PRC examples), `compute_bigG_from_g` folds it into the
two-event map by linear interpolation, and `simulate_rtm_plot_g` chains the
two for direct comparison with the abstract maps above.

## What to look for

For every two-event map plot, find intersections with the diagonal and then
compare the nearby slope with one. Iterate the two-oscillator cases to confirm
which intersections attract trajectories. Comparing the symmetric and asymmetric
abstract cases shows that symmetry of a phase relation is a model property, not
a general guarantee. Use `simulate_rtm_plot_g` to connect the map to a neuron model.

## Suggested order

1. Run `g_pulse_1` through `g_pulse_3` with the shared plotting helpers.
2. Compare `g_pulse_4`, `g_pulse_5`, and `simulate_f_tilde`.
3. Run both `simulate_two_pulse_coupled_osc` examples, then
   `simulate_rtm_plot_g`.

## Prerequisites and related chapters

Chapter 25 defines phase responses and interaction functions. Chapter 27 adds
delays to the pulse-map construction, and Chapter 28 develops a continuous
weak-coupling version of phase-difference dynamics.

## Running the examples

Open [`chapter26.ipynb`](chapter26.ipynb) in Jupyter, or via the Colab badge
at the top of the notebook, and run all cells top to bottom. The abstract maps
are quick NumPy/Matplotlib plots; `simulate_rtm_plot_g` integrates many RTM
phases and takes longer.
