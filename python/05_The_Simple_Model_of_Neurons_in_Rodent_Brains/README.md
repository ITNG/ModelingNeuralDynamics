# Simple models of neurons in rodent brains

## Overview

This chapter compares the RTM excitatory-neuron model with Wang--Buzsáki (WB)
and Erisir inhibitory-neuron models. All are conductance-based, but their
maximal conductances, reversal potentials, and gating kinetics produce
distinct spike waveforms and recovery speeds.

## Core ideas

The models share sodium, potassium, and leak currents, while their rate laws
encode cell-type-specific kinetics. In these implementations, sodium
activation \(m\) is set instantaneously to its steady-state value; \(h\) and
\(n\) remain dynamic. The two Erisir traces differ in the power used for the
potassium gate, making that modeling choice directly visible.

## Essential model

The common current-balance form is

\[
C\frac{dV}{dt}=I_{\mathrm{ext}}-g_{\mathrm{Na}}m_\infty(V)^3h(V-E_{\mathrm{Na}})
-g_{\mathrm{K}}n^p(V-E_{\mathrm{K}})-g_{\mathrm{L}}(V-E_{\mathrm{L}}).
\]

Here \(V\) is membrane voltage, \(t\) is time, \(C\) is capacitance,
\(I_{\mathrm{ext}}\) is applied current, \(g_{\mathrm{Na}}\),
\(g_{\mathrm{K}}\), and \(g_{\mathrm{L}}\) are maximal conductances,
\(E_{\mathrm{Na}}\), \(E_{\mathrm{K}}\), and \(E_{\mathrm{L}}\) are reversal
potentials, \(m_\infty\) is the instantaneous sodium-activation value, \(h\)
is sodium inactivation, \(n\) is potassium activation, and \(p\) is the
potassium-gate exponent selected by the model. Dynamic gates obey
\(dx/dt=\alpha_x(V)(1-x)-\beta_x(V)x\), where \(x\) is a gate and
\(\alpha_x\), \(\beta_x\) are its voltage-dependent rates.

## Code examples

- [`RTM_VOLTAGE_TRACE`](RTM_VOLTAGE_TRACE/) runs the RTM model and saves
  `fig_5_2.png`.
- [`WB_VOLTAGE_TRACE`](WB_VOLTAGE_TRACE/) runs the WB trace in `main.py`,
  saving `fig_5_3.png`; the alternate `main_wb_1996.py` entry point saves
  `fig_5_3_1996.png`.
- [`ERISIR_VOLTAGE_TRACE`](ERISIR_VOLTAGE_TRACE/) uses the Erisir model with
  an \(n^2\) potassium term and saves `fig_5_4.png`.
- [`ERISIR_VOLTAGE_TRACE_2`](ERISIR_VOLTAGE_TRACE_2/) changes that potassium
  term to \(n^4\) and saves `fig_5_4_2.png`.
- [`THREE_MODELS_GATING_VARIABLES`](THREE_MODELS_GATING_VARIABLES/) contains
  `RTM_GATING_VARIABLES/main.py`, which plots RTM gate steady states and time
  constants as `fig_5_1.png`.

## What to look for

Compare voltage traces before treating a label such as “excitatory” or
“inhibitory” as a complete explanation. In the gating figure, compare both
the equilibrium curves and their time constants; rate laws, not merely the
names of the gates, distinguish the models.

## Suggested order

1. Begin with `THREE_MODELS_GATING_VARIABLES` and the RTM rate functions.
2. Run `RTM_VOLTAGE_TRACE`, then `WB_VOLTAGE_TRACE`.
3. Compare the two Erisir variants and locate the changed power of \(n\).

## Prerequisites and related chapters

Chapters 01 and 03 introduce conductance-based HH-style models and gate
kinetics. Chapter 09 extends the RTM model with slow adaptation currents.

## Running the examples

Run each `main.py` from its example directory with `python main.py`. The
alternate WB implementation can be run with `python main_wb_1996.py` from
`WB_VOLTAGE_TRACE`. The nested gating example is run from
`THREE_MODELS_GATING_VARIABLES/RTM_GATING_VARIABLES` with `python main.py`.
All examples require SciPy, NumPy, and Matplotlib/Pylab and save PNGs in their
current working directory.
