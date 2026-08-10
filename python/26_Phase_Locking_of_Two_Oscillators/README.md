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

If \(g\) is a phase-resetting curve, the scripts first form the single-pulse
reset map \(f(\phi)=\phi+g(\phi)\). They then account for the other
oscillator's complementary phase with \(F(\phi)=f(1-\phi)\), and form the
two-event phase-difference map \(G(\phi)=F(F(\phi))\). Locking satisfies

\[
G(\phi_*)=\phi_*.
\]

For a one-dimensional map, the fixed point is locally stable when
\(|G'(\phi_*)|<1\). The scripts plot the two-event map and identity line so those
intersections can be read directly.

## Code examples

- [`ABSTRACT_PULSE_COUPLING_1`](ABSTRACT_PULSE_COUPLING_1/) plots an asymmetric
  positive abstract PRC and its phase map.
- [`ABSTRACT_PULSE_COUPLING_2`](ABSTRACT_PULSE_COUPLING_2/) changes the
  abstract response to show a different fixed-point arrangement.
- [`ABSTRACT_PULSE_COUPLING_3`](ABSTRACT_PULSE_COUPLING_3/) uses a PRC that is
  antisymmetric about phase one-half and gives its pulse-coupling map.
- [`ABSTRACT_PULSE_COUPLING_4`](ABSTRACT_PULSE_COUPLING_4/) uses an asymmetric
  PRC and maps its shifted locking geometry.
- [`ABSTRACT_PULSE_COUPLING_5`](ABSTRACT_PULSE_COUPLING_5/) supplies a further
  asymmetric response example.
- [`F_TILDE`](F_TILDE/) plots the transformed map used to locate phase-map
  fixed points.
- [`TWO_PULSE_COUPLED_OSC`](TWO_PULSE_COUPLED_OSC/) iterates two abstract
  pulse-coupled oscillators from selected phases.
- [`TWO_PULSE_COUPLED_OSC_2`](TWO_PULSE_COUPLED_OSC_2/) repeats the two-cell
  iteration for a second pulse-response choice.
- [`RTM_PLOT_G`](RTM_PLOT_G/) computes the RTM phase-resetting curve used for
  comparison with the abstract maps.

## What to look for

For every two-event map plot, find intersections with the diagonal and then
compare the nearby slope with one. Iterate the two-oscillator cases to confirm
which intersections attract trajectories. Comparing the symmetric and asymmetric
abstract cases shows that symmetry of a phase relation is a model property, not
a general guarantee. Use `RTM_PLOT_G` to connect the map to a neuron model.

## Suggested order

1. Run `ABSTRACT_PULSE_COUPLING_1` through `ABSTRACT_PULSE_COUPLING_3`.
2. Compare `ABSTRACT_PULSE_COUPLING_4`, `ABSTRACT_PULSE_COUPLING_5`, and
   `F_TILDE`.
3. Run `TWO_PULSE_COUPLED_OSC`, `TWO_PULSE_COUPLED_OSC_2`, and `RTM_PLOT_G`.

## Prerequisites and related chapters

Chapter 25 defines phase responses and interaction functions. Chapter 27 adds
delays to the pulse-map construction, and Chapter 28 develops a continuous
weak-coupling version of phase-difference dynamics.

## Running the examples

Run `python main.py` from an immediate example directory. The abstract maps are
quick NumPy/Matplotlib plots; `RTM_PLOT_G` integrates many RTM phases and takes
longer.
