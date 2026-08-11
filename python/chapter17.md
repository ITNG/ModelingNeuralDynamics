# Frequency-current curves

## Overview

These scripts compute firing rate as a function of applied current across
LIF, theta, HH, RTM, Erisir, reduced HH, INaP-I$_K$, and self-exciting theta
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

$$
T=\tau_m\log\!\left(\frac{\tau_m I}{\tau_m I-1}\right),\qquad f=1000/T.
$$

Here $I$ is normalized applied current, $\tau_m$ is membrane time
constant, $T$ is in milliseconds, and $f$ is in hertz. The conductance
models numerically count threshold crossings after transients.

## Code examples

All ten primary examples (plus three legacy full-model variants) now live
in one notebook, [`chapter17.ipynb`](chapter17.ipynb): `simulate_lif_f_i_curve`
and `simulate_theta_f_i_curve` plot the analytic LIF and theta-neuron f-I
curves; `simulate_setn_f_i` computes the self-exciting theta f-I curve;
`simulate_hh_reduced_f_i_curve`, `simulate_erisir_f_i_curve`, and
`simulate_wb_f_i_curve` compute forward/backward reduced-HH, Erisir, and
Wang-Buzsaki rates; `simulate_wb_f_i_curve_at_onset` magnifies WB onset
behavior; `simulate_rtm_with_m_current_f_i` measures the adapting RTM
model's f-I curve; `simulate_inapik_f_i_curve` computes INaP-I$_K$
forward/backward branches; `simulate_inapik_saddle_cycle_distance` plots
the distance between an INaP-I$_K$ saddle and cycle. A "Legacy Full-Model
F-I Curves" section preserves the original (untested) full HH and RTM RK4
scans as `simulate_hh_f_i_curve_legacy`, `simulate_rtm_f_i_curve_legacy`,
and `simulate_rtm_f_i_curve_at_onset_legacy`.

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
the INaP-I$_K$ and self-exciting theta systems; Chapter 18 explains the
bistability visible in some forward/backward curves.

## Running the examples

Open [`chapter17.ipynb`](chapter17.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom. Most
scans are a forward+backward sweep over dozens of currents and take tens
of seconds to a few minutes each; the notebook notes the slower ones
inline.
