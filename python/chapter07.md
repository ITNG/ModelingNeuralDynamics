# Linear integrate-and-fire neurons

## Overview

The linear integrate-and-fire (LIF) model keeps the passive accumulation of a
membrane but replaces an action potential's ionic mechanism with a threshold
and reset rule. The examples also measure the effective subthreshold behavior
of the HH model before comparing the two descriptions.

## Core ideas

Below threshold, voltage decays toward a driven equilibrium with membrane time
constant $\tau_m$. Crossing the chosen threshold is an event: the script
records a reset instead of resolving spike-generating sodium and potassium
currents. This makes LIF economical, but it cannot reproduce the HH gates or
their biophysical refractory dynamics without additional rules.

## Essential model

The subthreshold equation and event rule are

$$
\frac{dV}{dt}=-\frac{V}{\tau_m}+I,\qquad
V\ge V_{\mathrm{th}}\Rightarrow V\leftarrow V_{\mathrm{reset}}.
$$

Here $V$ is the model voltage, $t$ is time, $\tau_m$ is the membrane
time constant, $I$ is constant input in the code's normalized units,
$V_{\mathrm{th}}$ is threshold, and $V_{\mathrm{reset}}$ is the voltage
assigned after a threshold crossing. For the HH comparison, the instantaneous
effective time constant is $\tau=C/(g_{\mathrm{K}}n^4+g_{\mathrm{Na}}m^3h+g_{\mathrm{L}})$,
where $C$ is capacitance, $g_{\mathrm{K}}$, $g_{\mathrm{Na}}$, and
$g_{\mathrm{L}}$ are conductances, and $m$, $n$, and $h$ are HH gates.

## Code examples

All five examples now live in one notebook, [`chapter07.ipynb`](chapter07.ipynb),
with two `simulate_*` functions covering them: `simulate_hh_subthreshold`
backs `LIF_NEURON_WITH_HH`, `SUBTHR_FOR_HH`, and `TAU_M_FOR_HH`'s figures
(HH voltage vs. the LIF linear approximation, the current decomposition,
and the effective time constant); `simulate_lif_neuron(tau_m, i)` backs
`LIF_VOLTAGE_TRACE` and `LIF_VOLTAGE_TRACE_2`. Each has an `ipywidgets`
slider to explore its parameters interactively.

## What to look for

In the LIF traces, the discontinuous reset is a declared convention, not a
numerically resolved membrane event. Compare it with the smooth HH trace in
`LIF_NEURON_WITH_HH`, then use `SUBTHR_FOR_HH` and `TAU_M_FOR_HH` to see why a
single fixed passive time constant is only an approximation.

## Suggested order

1. Run `LIF_VOLTAGE_TRACE` and identify threshold and reset.
2. Run `LIF_VOLTAGE_TRACE_2` to connect input strength with interspike timing.
3. Read the three HH comparison examples in the order listed above.

## Prerequisites and related chapters

Chapters 01--04 provide the HH currents, gates, and trajectories used for
comparison. Chapter 08 replaces the linear subthreshold term with a quadratic
one; Chapter 09 adds slow adaptation to a LIF-style model.

## Running the examples

Open [`chapter07.ipynb`](chapter07.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook. Run all cells top to bottom; each
section's static figure reproduces the book's plot, and the `interact(...)`
cell below it lets you adjust that example's parameters with sliders.
