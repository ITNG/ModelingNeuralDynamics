# Spike-timing-dependent plasticity (STDP)

## Overview

STDP changes synaptic weights according to relative pre/post spike timing. The examples begin with the Abbott--Song rule and an adapting RTM trace, then progress through three-cell PING and a PING network with evolving excitatory synapses.

## Core ideas

Weight-update sign and size depend on spike order and time difference. An adaptation variable changes excitability and can alter order. In three-cell PING, E-to-E coupling changes E-assembly lag and frequency; active STDP makes that timing and the weights coevolve.

## Essential model

For \(\Delta t=t_{\rm post}-t_{\rm pre}\), a representative rule is \(\Delta w=A_+e^{-\Delta t/\tau_+}\) for post-after-pre and \(\Delta w=-A_-e^{\Delta t/\tau_-}\) for the opposite order. The PING implementation bounds recurrent E-to-E weights while spiking continues.

## Code examples

- [`ABBOTT_SONG`](ABBOTT_SONG/) plots the Abbott--Song timing-dependent rule.
- [`PING_WITH_STDP`](PING_WITH_STDP/) simulates a PING network with evolving recurrent E-to-E weights.
- [`RTM_VOLTAGE_TRACE_WITH_A`](RTM_VOLTAGE_TRACE_WITH_A/) shows RTM voltage and its adaptation variable.
- [`THREE_CELL_PING_1`](THREE_CELL_PING_1/) is the baseline two-E, one-I PING progression and writes spike times.
- [`THREE_CELL_PING_2`](THREE_CELL_PING_2/) compares weak and strong reciprocal E-to-E coupling.
- [`THREE_CELL_PING_3`](THREE_CELL_PING_3/) sweeps reciprocal E-to-E strength and writes lag/frequency data.
- [`THREE_CELL_PING_4`](THREE_CELL_PING_4/) sweeps one-way E-to-E strength and writes lag/frequency data.
- [`THREE_CELL_PING_5`](THREE_CELL_PING_5/) integrates three-cell PING with bounded STDP E-to-E updates.

## What to look for

Read the Abbott--Song curve with signed pre/post timing before interpreting weight traces. In three-cell cases, compare E-cell lag and frequency as coupling changes. In STDP cases, ask whether weight evolution changes the spike order that produced it.

## Suggested order

1. Run `ABBOTT_SONG` and `RTM_VOLTAGE_TRACE_WITH_A`.
2. Run `THREE_CELL_PING_1` through `THREE_CELL_PING_4`.
3. Finish with `THREE_CELL_PING_5` and `PING_WITH_STDP`.

## Prerequisites and related chapters

Chapter 9 provides adaptation context, Chapter 20 chemical synapses, Chapter 30 PING, and Chapter 39 short-term plasticity.

## Running the examples

Run `python main.py` in an immediate example directory. Some three-cell examples write spike-time or sweep data for auxiliary plotting scripts. NumPy, SciPy, and Matplotlib are required.
