# Linear integrate-and-fire neurons

## Overview

The linear integrate-and-fire (LIF) model keeps the passive accumulation of a
membrane but replaces an action potential's ionic mechanism with a threshold
and reset rule. The examples also measure the effective subthreshold behavior
of the HH model before comparing the two descriptions.

## Core ideas

Below threshold, voltage decays toward a driven equilibrium with membrane time
constant \(\tau_m\). Crossing the chosen threshold is an event: the script
records a reset instead of resolving spike-generating sodium and potassium
currents. This makes LIF economical, but it cannot reproduce the HH gates or
their biophysical refractory dynamics without additional rules.

## Essential model

The subthreshold equation and event rule are

\[
\frac{dV}{dt}=-\frac{V}{\tau_m}+I,\qquad
V\ge V_{\mathrm{th}}\Rightarrow V\leftarrow V_{\mathrm{reset}}.
\]

Here \(V\) is the model voltage, \(t\) is time, \(\tau_m\) is the membrane
time constant, \(I\) is constant input in the code's normalized units,
\(V_{\mathrm{th}}\) is threshold, and \(V_{\mathrm{reset}}\) is the voltage
assigned after a threshold crossing. For the HH comparison, the instantaneous
effective time constant is \(\tau=C/(g_{\mathrm{K}}n^4+g_{\mathrm{Na}}m^3h+g_{\mathrm{L}})\),
where \(C\) is capacitance, \(g_{\mathrm{K}}\), \(g_{\mathrm{Na}}\), and
\(g_{\mathrm{L}}\) are conductances, and \(m\), \(n\), and \(h\) are HH gates.

## Code examples

- [`SUBTHR_FOR_HH`](SUBTHR_FOR_HH/) calculates the HH sodium, potassium,
  leak, and total subthreshold currents and saves `fig_7_3.png`.
- [`TAU_M_FOR_HH`](TAU_M_FOR_HH/) computes the HH effective membrane time
  constant during a trajectory and saves `fig_7_2.png`.
- [`LIF_VOLTAGE_TRACE`](LIF_VOLTAGE_TRACE/) uses RK4 integration with a
  threshold of 1 and reset of 0, saving `fig_7_4.png`.
- [`LIF_VOLTAGE_TRACE_2`](LIF_VOLTAGE_TRACE_2/) adjusts the drive to target
  its timing and saves `fig_7_5.png`.
- [`LIF_NEURON_WITH_HH`](LIF_NEURON_WITH_HH/) overlays reset-like LIF ramps
  with a full HH voltage trace and saves `fig_7_1.png`.

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

Run `python main.py` from any example directory. `LIF_VOLTAGE_TRACE_2` opens
its Matplotlib window after writing `fig_7_5.png`; the remaining examples save
their PNGs without an active `show()` call. SciPy is required for the HH
comparison examples, while the standalone LIF traces use NumPy and Matplotlib.
