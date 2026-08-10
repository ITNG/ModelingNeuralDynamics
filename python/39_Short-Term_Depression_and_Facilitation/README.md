# Short-term depression and facilitation

## Overview

Synapses can change strength through a pulse train even with fixed anatomy.
These examples first introduce a smooth spike/pulse approximation motivating
synaptic gating, then show depressing and facilitating synapses in RTM and WB
neurons.

## Core ideas

Depression depletes available resources after release and recovers between events. Facilitation transiently increases utilization, so a pulse sequence can initially grow. Their balance changes effective synaptic conductance, postsynaptic timing, and whether early or late pulses evoke spikes.

## Essential model

With resources \(x\) and utilization \(u\), release is proportional to
\(ux\); depression lowers \(x\) after events and facilitation raises \(u\).
The resulting synaptic gate enters the conductance current
\(gs(E_{\rm rev}-v)\). The introductory smooth proxy is
\(\gamma=1+\tanh(v/10)\), which approximates a spike-shaped pulse for a
voltage trace.

## Code examples

- [`PULSES`](PULSES/) plots an RTM voltage trace and its smooth
  \(\gamma=1+\tanh(v/10)\) spike/pulse approximation, including its
  per-spike area calculation.
- [`RTM_WITH_DEPRESSING_AND_FACILITATING_S`](RTM_WITH_DEPRESSING_AND_FACILITATING_S/) compares both mechanisms in RTM.
- [`RTM_WITH_DEPRESSING_S`](RTM_WITH_DEPRESSING_S/) shows depression alone in RTM.
- [`WB_WITH_DEPRESSING_S`](WB_WITH_DEPRESSING_S/) gives the depressing-synapse case for WB.

## What to look for

In `PULSES`, compare the voltage crossing with the smooth \(\gamma\) pulse;
it is not a resource/utilization simulation. In the remaining examples, track
synaptic strength over events: depression should weaken later events until
recovery, while facilitation can enhance them transiently. Compare RTM and WB
to separate synaptic from intrinsic-neuron effects.

## Suggested order

1. Run `PULSES` to see the smooth spike/pulse approximation.
2. Compare the depressing RTM and WB examples.
3. Run the RTM depression-plus-facilitation comparison.

## Prerequisites and related chapters

Chapter 20 supplies the chemical-synapse gate, Chapter 5 RTM and WB cells, and Chapter 40 spike-timing-dependent plasticity.

## Running the examples

Run `python main.py` from an immediate example directory. NumPy, SciPy, and Matplotlib are required.
