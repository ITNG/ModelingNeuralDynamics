# Quadratic integrate-and-fire and theta neurons

## Overview

Quadratic integrate-and-fire (QIF) neurons replace LIF's linear drift with a
quadratic voltage flow. Their apparent divergence is converted into a regular
spike cycle by a reset convention, or equivalently represented as smooth
motion around the circle by the theta transformation.

## Core ideas

A QIF trajectory can reach infinity in finite time; this is the model's spike
event, not a claim that a physical membrane voltage is infinite. Resetting
from the upper to lower branch makes a repeated trace. Under a voltage-to-phase
change of variables, the two infinities meet at one point on a circle, so the
same firing process becomes a continuous phase flow.

## Essential model

The QIF examples use

$$
\frac{dV}{dt}=-\frac{V}{\tau_m}(1-V)+I,
\qquad V=1-\cos\theta.
$$

Here $V$ is normalized voltage, $t$ is time, $\tau_m$ is the membrane
time constant, $I$ is constant input, and $\theta$ is circular phase.
The QIF trace resets after passing the finite numerical threshold $V=1$,
while the analytic infinite-threshold construction connects the positive and
negative infinite branches. In the circle diagram, the transition occurs at
$I=1/(4\tau_m)$, where $I$ and $\tau_m$ have the meanings above.

## Code examples

All four examples now live in one notebook, [`chapter08.ipynb`](chapter08.ipynb):
`simulate_qif_voltage_trace` integrates a QIF trace with a midpoint step and
reset convention; `simulate_qif_infinite_threshold` constructs the analytic
finite-time blow-up branches and their continuation; `simulate_theta_firing`
integrates phase $\theta$ and plots $1-\cos\theta$; `plot_three_circles`
draws phase flows below, at, and above the firing transition on three
circles (a schematic, no simulation). Each simulated example has an
`ipywidgets` slider to explore its parameters interactively.

## What to look for

The dashed segments in `QIF_VOLTAGE_TRACE` mark the reset convention. Compare
that with `QIF_INFINITE_THRESHOLD`, where the finite-time divergence is shown
explicitly. `THETA_FIRING` has no reset jump because phase wraps naturally;
`THREE_CIRCLES` explains why the three input regimes have different geometry.

## Suggested order

1. Run `QIF_VOLTAGE_TRACE` to see the practical event rule.
2. Run `QIF_INFINITE_THRESHOLD` to see what that rule abbreviates.
3. Run `THETA_FIRING`, then use `THREE_CIRCLES` to interpret the phase flow.

## Prerequisites and related chapters

Chapter 07 supplies the linear integrate-and-fire baseline. Familiarity with
one-dimensional ODEs, trigonometric functions, and phase portraits is useful.
Later network examples use theta-neuron phase descriptions.

## Running the examples

Open [`chapter08.ipynb`](chapter08.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook. Run all cells top to bottom; each
section's static figure reproduces the book's plot, and the `interact(...)`
cell below it lets you adjust that example's parameters with sliders.
