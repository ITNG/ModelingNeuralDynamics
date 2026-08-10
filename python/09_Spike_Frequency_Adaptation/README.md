# Spike-frequency adaptation

## Overview

These examples add slow negative feedback to spiking models. They cover a
voltage-dependent M-current, a calcium-dependent afterhyperpolarization (AHP)
current, a reset-based adaptation variable, and maps that describe how the
slow variable changes from one spike to the next.

## Core ideas

Adaptation accumulates during activity and reduces subsequent excitability,
which can lengthen interspike intervals. The M-current uses a slow gate,
whereas the AHP current uses a calcium-like state. Resting variants set the
external drive to zero, showing the models' zero-drive resting trajectories
alongside sustained spiking regimes. A spike-to-spike map condenses continuous
evolution into a one-dimensional update.

## Essential model

For the adaptation LIF example, the equations and spike update are

\[
\frac{dV}{dt}=-\frac{V}{\tau_m}+I-wV,\qquad
\frac{dw}{dt}=-\frac{w}{\tau_a},\qquad
V\ge1\Rightarrow (V,w)\leftarrow(0,w+\Delta).
\]

Here \(V\) is normalized voltage, \(t\) is time, \(\tau_m\) is membrane
time constant, \(I\) is applied input, \(w\) is the nonnegative adaptation
state, \(\tau_a\) is its decay time constant, and \(\Delta\) is its increment
at a spike. In the RTM M-current model, the added current is
\(I_M=g_Mw(V-E_K)\), where \(g_M\) is maximal M conductance and \(E_K\) is
the potassium reversal potential. In the AHP model it is
\(I_{AHP}=g_{AHP}c(V-E_K)\), where \(g_{AHP}\) is maximal AHP conductance and
\(c\) is the calcium-like state.

## Code examples

- [`M_CURRENT`](M_CURRENT/) plots the M-gate steady state and time constant,
  saving `fig.png`.
- [`CALCIUM_RISE`](CALCIUM_RISE/) plots the voltage-dependent calcium target
  used by the AHP examples, saving `fig.png`.
- [`RTM_M`](RTM_M/) adds an M-current to a driven RTM neuron and saves voltage
  and gate traces as `fig.png`.
- [`RTM_M_RESTING`](RTM_M_RESTING/) uses zero external input for the M-current
  model and saves the resting trajectory as `fig.png`.
- [`RTM_AHP`](RTM_AHP/) adds calcium-dependent AHP feedback to a driven RTM
  neuron and saves voltage and calcium traces as `fig.png`.
- [`RTM_AHP_RESTING`](RTM_AHP_RESTING/) runs the AHP model at zero input and
  saves its resting trajectory as `fig.png`.
- [`LIF_ADAPT`](LIF_ADAPT/) integrates the reset-based adaptation model and
  saves voltage and adaptation state as `fig.png`.
- [`ADAPTATION_MAP`](ADAPTATION_MAP/) computes the spike-to-spike map
  \(\phi(z)\), its diagonal, and a marked fixed point, saving `fig.png`.
- [`V_V_TILDE`](V_V_TILDE/) compares two subthreshold voltages with different
  initial adaptation amplitudes and saves `fig.png`.

## What to look for

In `RTM_M` and `RTM_AHP`, watch the slow state rise while spiking continues;
the resting versions show what remains when drive is removed. `LIF_ADAPT`
makes the discrete spike increment explicit. In `ADAPTATION_MAP`, intersections
of \(\phi(z)\) with the diagonal identify fixed spike-to-spike adaptation
levels, where \(z\) is the adaptation value just after a spike and
\(\phi(z)\) is its value after the next spike.

## Suggested order

1. Run `M_CURRENT` and `CALCIUM_RISE` to inspect the slow feedback laws.
2. Compare each driven RTM model with its `_RESTING` counterpart.
3. Run `LIF_ADAPT`, `V_V_TILDE`, and finally `ADAPTATION_MAP` to connect
   continuous slow dynamics with a per-spike description.

## Prerequisites and related chapters

Chapter 05 introduces the RTM neuron used here, and Chapter 07 supplies the
LIF reset convention. The QIF and theta models in Chapter 08 provide other
ways to represent spike events and phase.

## Running the examples

Run `python main.py` from each immediate example directory. The RTM examples
require SciPy, NumPy, and Matplotlib/Pylab; the remaining examples use NumPy
and Matplotlib. Every entry point writes `fig.png` in its working directory.
