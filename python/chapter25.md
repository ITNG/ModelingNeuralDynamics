# Phase response curves (PRCs)

## Overview

A phase response curve reports how much a single perturbation shifts a
neuron's firing phase: deliver it at phase $\varphi$ (fraction of the cycle
since the last spike) and measure the resulting shift $g(\varphi)$. The
examples move from the theta neuron's closed-form PRC to numerically
measured PRCs for RTM, WB, HH, and Erisir neurons, under both synaptic-pulse
and instantaneous-voltage-kick perturbations.

## Core ideas

For conductance-based neurons, $g$ is measured by splay-initializing a
population of independent neurons across every phase, delivering the same
perturbation to all of them at time 0, and recording each neuron's next
spike time $t_\ast$: $g(\varphi)=1-\varphi-t_\ast/T$, with $T$ the unperturbed
period. Weak-perturbation PRCs (small synaptic conductance or small voltage
kick) approach a linear-response limit -- normalizing by the perturbation
strength collapses curves at different strengths onto one shape.

## Essential model

Each PRC example reuses one of the conductance-based neuron models from
earlier chapters (RTM, WB, HH, Erisir), perturbed either by a synaptic pulse

$$
I_{\rm syn}=g_{\rm syn}\,s\,(v_{\rm rev}-v),\qquad s(0)=0,\ q(0)=1,
$$

or by an instantaneous voltage kick $v\to v+\Delta v$ at $t=0$. The theta
neuron instead has a closed-form return map
$f(\varphi)=\frac1\pi\arctan\!\big(\tan(\pi\varphi-\tfrac\pi2)+\delta v/\sqrt{\tau_mI-\tfrac14}\big)+\tfrac12$
and PRC $g=f-\varphi$.

## Code examples

All thirteen examples now live in one notebook, [`chapter25.ipynb`](chapter25.ipynb):
`theta_f`/`theta_prc`/`theta_prc_short_weak` give the theta neuron's closed-form
return map, PRC, and linear-response PRC. `rtm_init` (shared with Chapter 24)
finds the single-cell RTM limit cycle for splay initialization.
`simulate_synaptic_pulse_prc` computes the RTM synaptic-pulse PRC, reused
directly for the interaction-function panel ($f=\varphi+g$), the three-weak-pulses
comparison (`simulate_three_weak_prcs`), and the weak-pulse linear-response
check (`simulate_weak_prc_comparison`). `simulate_voltage_kick_prc` is the
voltage-kick analogue, reused for the strong-kick and weak-kick-comparison
(`simulate_short_weak_prc_comparison`) examples. `simulate_phase_shift` shows
the shift directly as two overlaid voltage traces (`simulate_baseline` and
`simulate_perturbed`) instead of a PRC curve. `simulate_wb_prc` (shared WB
gating and `wb_init`) covers both the inhibitory-pulse WB example and the WB
panel of the model-comparison figure; `simulate_hh_prc` and `simulate_erisir_prc`
provide the other two panels of that comparison.

## What to look for

Compare the theta neuron's smooth, everywhere-positive PRC to the
conductance-based neurons' PRCs, which can change sign (a pulse can delay or
advance firing depending on phase) or even go negative across most of the
cycle for inhibition. Check that halving the perturbation strength collapses
the normalized PRC onto the same curve (linear response), and that the
phase-shift trace directly shows the same effect the PRC curve predicts.

## Suggested order

1. Start with the theta-neuron trio (`theta_f`, `theta_prc`,
   `theta_prc_short_weak`) for the closed-form baseline.
2. Move to `simulate_synaptic_pulse_prc` and the interaction-function panel,
   then the three-weak-pulses and weak-pulse linear-response checks.
3. Compare with `simulate_voltage_kick_prc` and its weak-kick check, then
   `simulate_phase_shift` for the same effect shown as raw traces.
4. Finish with the WB inhibitory-pulse example and the WB/HH/Erisir
   model-comparison panel.

## Prerequisites and related chapters

Chapter 20 develops the synaptic gate kinetics ($\tau_r,\tau_d,\tau_{dq}$)
used throughout. Chapter 24's RTM limit-cycle finder is reused directly here.
Chapters 26-29 build on PRCs to analyze phase locking and synchrony stability
in coupled-oscillator networks.

## Running the examples

Open [`chapter25.ipynb`](chapter25.ipynb) in Jupyter, or via the Colab badge
at the top of the notebook, and run all cells top to bottom. Several sliders
(`interact`) let you adjust perturbation strength or timing and watch the
PRC or trace update.
