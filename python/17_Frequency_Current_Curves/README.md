# Frequency-current curves

## Overview

These scripts compute firing rate as a function of applied current across
LIF, theta, HH, RTM, Erisir, reduced HH, INaP-I\(_K\), and self-exciting theta
models. Forward and backward scans reveal onset frequency, continuous or
discontinuous branches, and bistability.

## Core ideas

An f--I curve maps steady firing rate to input current. Type-1 onset rises
continuously from zero frequency, while type-2 onset starts at a nonzero
frequency. Carrying the terminal state from one current to the next gives a
forward branch; reversing the scan can reveal a different branch when resting
and spiking attractors coexist.

## Essential model

For the LIF example, the interspike interval and firing rate are

\[
T=\tau_m\log\!\left(\frac{\tau_m I}{\tau_m I-1}\right),\qquad f=1000/T.
\]

Here \(I\) is normalized applied current, \(\tau_m\) is membrane time
constant, \(T\) is in milliseconds, and \(f\) is in hertz. The conductance
models numerically count threshold crossings after transients.

## Code examples

- [`LIF_F_I_CURVE`](LIF_F_I_CURVE/) plots the analytic LIF f--I curve in
  `fig.png`.
- [`THETA_F_I_CURVE`](THETA_F_I_CURVE/) plots the theta-neuron f--I relation in
  `fig.png`.
- [`SETN_F_I`](SETN_F_I/) computes the self-exciting theta f--I curve in
  `fig.png`.
- [`HH_F_I_CURVE`](HH_F_I_CURVE/) runs forward and backward HH scans with
  `HH_F_I_CURVE.py`, saves data files, and writes `fig_17_1.png`.
- [`HH_REDUCED_F_I_CURVE`](HH_REDUCED_F_I_CURVE/) computes forward and backward
  reduced-HH rates in `fig.png`.
- [`RTM_F_I_CURVE`](RTM_F_I_CURVE/) performs the RTM scan with
  `RTM_F_I_CURVE.py` and saves `fig_17_1.png`.
- [`RTM_F_I_CURVE_AT_ONSET`](RTM_F_I_CURVE_AT_ONSET/) resolves RTM onset with
  `RTM_F_I_CURVE_AT_ONSET.py` and saves `fig_17_5.png`.
- [`RTM_WITH_M_CURRENT_F_I`](RTM_WITH_M_CURRENT_F_I/) measures the adapting
  RTM model's f--I curve in `fig.png`.
- [`ERISIR_F_I_CURVE`](ERISIR_F_I_CURVE/) computes Erisir forward and backward
  branches in `fig.png`.
- [`INAPIK_F_I_CURVE`](INAPIK_F_I_CURVE/) computes INaP-I\(_K\) forward and
  backward branches in `fig.png`.
- [`INAPIK_SADDLE_CYCLE_DISTANCE`](INAPIK_SADDLE_CYCLE_DISTANCE/) plots the
  distance between an INaP-I\(_K\) saddle and cycle in `fig.png`.
- [`WB_F_I_CURVE`](WB_F_I_CURVE/) computes the Wang--Buzsaki f--I curve in
  `fig.png`.
- [`WB_F_I_CURVE_AT_ONSET`](WB_F_I_CURVE_AT_ONSET/) magnifies Wang--Buzsaki
  onset behavior in `fig.png`.

## What to look for

Compare `LIF_F_I_CURVE` and `THETA_F_I_CURVE` for continuous zero-frequency
onset. Use each paired forward/backward scan to find hysteresis. The two
`_AT_ONSET` examples use finer current ranges, and
`INAPIK_SADDLE_CYCLE_DISTANCE` links a geometric distance to the transition.

## Suggested order

1. Run `LIF_F_I_CURVE`, `THETA_F_I_CURVE`, and `SETN_F_I`.
2. Compare `HH_REDUCED_F_I_CURVE`, `ERISIR_F_I_CURVE`, and `INAPIK_F_I_CURVE`.
3. Use the RTM, HH, and WB onset scripts for detailed model comparisons.

## Prerequisites and related chapters

Chapters 12 and 14 distinguish type-1 and type-2 onset. Chapter 16 provides
the INaP-I\(_K\) and self-exciting theta systems; Chapter 18 explains the
bistability visible in some forward/backward curves.

## Running the examples

Run `python main.py` in directories that contain `main.py`. For the legacy HH
and RTM directories, run `python HH_F_I_CURVE.py`, `python RTM_F_I_CURVE.py`,
or `python RTM_F_I_CURVE_AT_ONSET.py` as appropriate; their `fff.py` scripts
replot saved scan data. These long scans use NumPy, SciPy, and Matplotlib.
