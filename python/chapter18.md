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

All sixteen examples now live in one notebook, [`chapter18.ipynb`](chapter18.ipynb):
`simulate_hh_bistable`/`_gates`/`_limited_n` compare resting and firing HH
trajectories, plot HH gates and their targets, and isolate the potassium
gate's role; `simulate_erisir_bistable`/`_gates`/`_limited_h` do the same
for the Erisir model; `simulate_h_current` and `simulate_modified_tau_r`
plot the h-gate steady state and (modified) time constant;
`simulate_rtm_voltage_trace_with_i_h` simulates a rebound-firing RTM
trace; `simulate_rtm_f_i_curve_with_i_h` computes h-current RTM
forward/backward f-I branches; `simulate_rtm_with_i_h_bistable`/`_gates`/
`simulate_rtm_with_i_h_limited_r` do the h-current RTM bistability family;
`simulate_rtm_with_i_m_bistable`/`_gates`/`simulate_rtm_with_i_m_limited_w`
do the M-current RTM family.

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

Open [`chapter18.ipynb`](chapter18.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom. The
`simulate_rtm_f_i_curve_with_i_h` cell is very slow (several minutes to
tens of minutes) -- the notebook notes this inline.
