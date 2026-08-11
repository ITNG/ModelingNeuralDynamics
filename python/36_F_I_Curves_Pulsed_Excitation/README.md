# f--I curves with pulsed excitation

## Overview

This chapter compares firing-rate curves driven by periodic excitation with steady-current f--I curves. Square pulses provide an idealized input, and RTM simulations show how pulse timing and amplitude select the spikes that determine rate.

## Core ideas

A pulsed f--I curve counts response opportunities per pulse period. It can have steps or saturation because a cell may fire once, skip, or fire more than once per pulse. The square-pulse model isolates this temporal effect before conductance-based RTM dynamics are added.

## Essential model

For pulse period $T$, the observed rate is $f=N_{\rm spikes}/T_{\rm observation}$. Unlike steady injected current, drive is nonzero only during a pulse interval; sweeping pulse amplitude produces the pulsed f--I curve.

## Code examples

- [`IDEALIZED_F_I_CURVE`](IDEALIZED_F_I_CURVE/) plots an idealized pulsed f--I construction.
- [`RTM_F_I_CURVE_PULSED_EXCITATION`](RTM_F_I_CURVE_PULSED_EXCITATION/) computes the RTM pulsed-excitation f--I curve.
- [`RTM_F_I_CURVE_PULSED_EXCITATION_2`](RTM_F_I_CURVE_PULSED_EXCITATION_2/) compares a second RTM pulse regime.
- [`SQUARE_PULSES`](SQUARE_PULSES/) plots the idealized square-pulse forcing signal.

## What to look for

Read rate changes with the pulse schedule: a jump can be an extra spike per pulse and a plateau can be a locked one-spike response. Compare the idealized curve with both RTM outputs before attributing all shape to intrinsic f--I nonlinearity.

## Suggested order

1. Run `SQUARE_PULSES` and `IDEALIZED_F_I_CURVE`.
2. Run `RTM_F_I_CURVE_PULSED_EXCITATION`.
3. Compare it with `RTM_F_I_CURVE_PULSED_EXCITATION_2`.

## Prerequisites and related chapters

Chapter 17 provides steady-current f--I curves, Chapter 23 periodic excitation and entrainment, and Chapter 35 periodic inhibition.

## Running the examples

Run `python main.py` in an immediate example directory. The scripts require NumPy and Matplotlib; RTM sweeps can take longer.
