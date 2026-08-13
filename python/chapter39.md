# Short-Term Depression and Facilitation

## Overview

Synapses can change strength through a pulse train even with fixed anatomy.
The examples in [`chapter39.ipynb`](chapter39.ipynb) first introduce a
smooth spike/pulse approximation motivating synaptic gating, then show
depressing and facilitating synapses in RTM and WB neurons.

## Core ideas

Depression depletes available resources after release and recovers between
events. Facilitation transiently increases utilization, so a pulse
sequence can initially grow. Their balance changes effective synaptic
conductance, postsynaptic timing, and whether early or late pulses evoke
spikes.

## Essential model

With resources $x$ and utilization $u$, release is proportional to $ux$;
depression lowers $x$ after events and facilitation raises $u$. The
resulting synaptic gate enters the conductance current $gs(E_{\rm
rev}-v)$. The introductory smooth proxy is $\gamma=1+\tanh(v/10)$, which
approximates a spike-shaped pulse for a voltage trace.

## Code examples

All four examples live in one notebook, [`chapter39.ipynb`](chapter39.ipynb).

`simulate_pulses`/`plot_pulses` cover `PULSES`: an RTM voltage trace and
its smooth $\gamma=1+\tanh(v/10)$ spike/pulse approximation, including its
per-spike area calculation, with `i_ext` exposed through an `interact()`
slider.

`simulate_rtm_depressing_facilitating`/`plot_rtm_depressing_facilitating`
cover `RTM_WITH_DEPRESSING_AND_FACILITATING_S`, comparing depression and
facilitation together in an RTM cell, with `tau_facil` exposed through a
slider.

`simulate_rtm_depressing`/`plot_rtm_depressing` cover
`RTM_WITH_DEPRESSING_S`, depression alone in an RTM cell, with `U` and
`tau_rec` exposed through sliders.

`simulate_wb_with_depressing_s`/`plot_wb_with_depressing_s` cover
`WB_WITH_DEPRESSING_S`, the depressing-synapse case for a WB cell
(integrated with an explicit Heun/RK2 step rather than `odeint`), with `U`
and `C` exposed through sliders.

## What to look for

In `PULSES`, compare the voltage crossing with the smooth $\gamma$ pulse;
it is not a resource/utilization simulation. In the remaining examples,
track synaptic strength over events: depression should weaken later
events until recovery, while facilitation can enhance them transiently.
Compare RTM and WB to separate synaptic from intrinsic-neuron effects.

## Suggested order

1. Run `PULSES` to see the smooth spike/pulse approximation.
2. Compare the depressing RTM and WB examples.
3. Run the RTM depression-plus-facilitation comparison.

## Prerequisites and related chapters

Chapter 20 supplies the chemical-synapse gate, Chapter 5 RTM and WB cells,
and Chapter 40 spike-timing-dependent plasticity.

## Running the examples

Open [`chapter39.ipynb`](chapter39.ipynb) and run the cells top to bottom.
NumPy, SciPy, Matplotlib, and ipywidgets are required.
