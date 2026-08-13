# f--I curves with pulsed excitation

## Overview

This chapter compares firing-rate curves driven by periodic excitation with
steady-current f--I curves. Square pulses provide an idealized input, and RTM
simulations show how pulse timing and amplitude select the spikes that
determine rate.

## Core ideas

A pulsed f--I curve counts response opportunities per pulse period. It can
have steps or saturation because a cell may fire once, skip, or fire more
than once per pulse. The square-pulse model isolates this temporal effect
before conductance-based RTM dynamics are added.

## Essential model

For pulse period $T$, the observed rate is $f=N_{\rm spikes}/T_{\rm
observation}$. Unlike steady injected current, drive is nonzero only during
a pulse interval; sweeping pulse amplitude produces the pulsed f--I curve.

## Code examples

All four examples now live in one notebook,
[`chapter36.ipynb`](chapter36.ipynb): `plot_idealized_f_i_curve` draws the
schematic type-1-onset $f=\sqrt{I-I_c}$ construction; `plot_square_pulses`
draws the idealized square-pulse forcing signal. The two RTM examples share
a single simulation kernel (`_rtm_step`, `_rtm_f_i_curve_constant_python`,
`_rtm_f_i_curve_pulsed_python`, `rtm_pulse_shape`) driven through
`_compute_rtm_f_i_curves`; `compute_rtm_f_i_curves_pulsed_excitation`
(`g_l=0.1`) and `compute_rtm_f_i_curves_pulsed_excitation_2` (`g_l=0.2`)
each return `(f_vec_constant, f_vec_pulsed, i_ext_vec)`, plotted with
`plot_rtm_f_i_curves_pulsed_excitation`.

## What to look for

Read rate changes with the pulse schedule: a jump can be an extra spike per
pulse and a plateau can be a locked one-spike response. Compare the
idealized curve with both RTM outputs before attributing all shape to
intrinsic f--I nonlinearity.

## Suggested order

1. Run the "Square Pulses" and "Idealized F-I Curve" sections.
2. Run "RTM F-I Curve, Pulsed Excitation".
3. Compare it with "RTM F-I Curve, Pulsed Excitation 2".

## Prerequisites and related chapters

Chapter 17 provides steady-current f--I curves, Chapter 23 periodic
excitation and entrainment, and Chapter 35 periodic inhibition.

## Running the examples

Open [`chapter36.ipynb`](chapter36.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom. The
schematic sections are instant; each RTM F-I sweep integrates 1000 ms
across 201 drive values and takes roughly ten seconds with `numba`.
