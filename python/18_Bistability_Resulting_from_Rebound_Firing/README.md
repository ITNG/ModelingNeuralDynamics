# Bistability resulting from rebound firing

## Overview

These examples examine coexistence of resting and spiking behavior, rebound
firing, and hysteresis caused by slow gates. They contrast HH and Erisir
models with RTM neurons carrying either h-current or M-current feedback.

## Core ideas

Bistability means the same drive supports more than one attracting behavior;
the initial condition selects rest or sustained firing. A slowly evolving
gate shifts effective excitability and can produce rebound after inhibition.
Forward and backward f--I sweeps expose hysteresis, while gate plots compare
actual gates to their voltage-dependent steady-state values.

## Essential model

The added h-current has the form

$$
I_h=g_hr(v_h-v),\qquad \dot r=\frac{r_\infty(v)-r}{\tau_r(v)},
$$

where $r$ is a slow activation gate, $g_h$ is maximal conductance, and
$v_h$ is its reversal potential. The M-current analog is
$I_M=g_Mw(v_K-v)$ with a slow potassium gate $w$.

## Code examples

- [`HH_BISTABLE`](HH_BISTABLE/) compares resting and firing HH trajectories in
  `fig.png`.
- [`HH_BISTABLE_GATES`](HH_BISTABLE_GATES/) plots HH gates and their targets in
  `fig.png`.
- [`HH_BISTABLE_LIMITED_N`](HH_BISTABLE_LIMITED_N/) limits the HH potassium
  gate and plots the voltage response in `fig.png`.
- [`ERISIR_BISTABLE`](ERISIR_BISTABLE/) compares rest and firing in the Erisir
  model in `fig.png`.
- [`ERISIR_BISTABLE_GATES`](ERISIR_BISTABLE_GATES/) plots Erisir gate dynamics
  in `fig.png`.
- [`ERISIR_BISTABLE_LIMITED_H`](ERISIR_BISTABLE_LIMITED_H/) limits Erisir
  sodium inactivation and saves `fig.png`.
- [`H_CURRENT`](H_CURRENT/) plots h-gate steady state and time constant in
  `fig.png`.
- [`PLOT_MODIFIED_TAU_R`](PLOT_MODIFIED_TAU_R/) compares original and modified
  h-gate time constants in `fig.png`.
- [`RTM_VOLTAGE_TRACE_WITH_I_H`](RTM_VOLTAGE_TRACE_WITH_I_H/) simulates an RTM
  voltage trace with h-current in `fig.png`.
- [`RTM_F_I_CURVE_WITH_I_H`](RTM_F_I_CURVE_WITH_I_H/) computes h-current RTM
  forward and backward f--I branches in `fig.png`.
- [`RTM_WITH_I_H_BISTABLE`](RTM_WITH_I_H_BISTABLE/) compares h-current RTM
  resting and firing trajectories in `fig.png`.
- [`RTM_WITH_I_H_BISTABLE_GATES`](RTM_WITH_I_H_BISTABLE_GATES/) shows its
  h-current and HH-like gates in `fig.png`.
- [`RTM_WITH_I_H_LIMITED_R`](RTM_WITH_I_H_LIMITED_R/) constrains the slow
  h-gate and plots voltage in `fig.png`.
- [`RTM_WITH_I_M_BISTABLE`](RTM_WITH_I_M_BISTABLE/) compares M-current RTM
  rest and firing in `fig.png`.
- [`RTM_WITH_I_M_BISTABLE_GATES`](RTM_WITH_I_M_BISTABLE_GATES/) plots its
  slow M gate with the fast gates in `fig.png`.
- [`RTM_WITH_I_M_LIMITED_W`](RTM_WITH_I_M_LIMITED_W/) constrains the M gate and
  saves the voltage response in `fig.png`.

## What to look for

Start each family with its `_BISTABLE` trace, then use `_GATES` to see which
slow variable differs from its instantaneous target. Compare h-current and
M-current RTM families: the former supplies depolarizing rebound feedback,
whereas the latter is a slow outward brake. The f--I curve identifies the
hysteretic interval.

## Suggested order

1. Run `H_CURRENT`, `PLOT_MODIFIED_TAU_R`, and `RTM_VOLTAGE_TRACE_WITH_I_H`.
2. Compare the three HH and three Erisir examples.
3. Compare the h-current and M-current RTM bistability families, then the f--I
   curve.

## Prerequisites and related chapters

Chapter 09 covers slow adaptation currents, Chapter 14 explains unstable
cycles and Chapter 17 introduces forward/backward f--I scans. Chapter 19 uses
slow-current feedback on an even longer time scale to generate bursting.

## Running the examples

Run `python main.py` from an immediate example directory. The scripts use
NumPy and Matplotlib and save `fig.png`; the f--I scan takes longer because it
integrates many currents to steady state.
