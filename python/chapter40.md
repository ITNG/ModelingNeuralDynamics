# Spike Timing-Dependent Plasticity (STDP)

## Overview

STDP changes synaptic weights according to the relative timing of pre-
and postsynaptic spikes. The examples in
[`chapter40.ipynb`](chapter40.ipynb) begin with the Abbott--Song timing
rule and an adapting RTM voltage trace, then progress through a
three-cell PING network (first with fixed E-to-E coupling, then with that
coupling evolving under STDP), and finish with a full PING network in
which every recurrent E-to-E synapse evolves under STDP.

## Core ideas

Weight-update sign and size depend on spike order and time difference: a
presynaptic spike shortly before a postsynaptic one potentiates the
synapse, and the reverse order depresses it. An adaptation-like variable
tracks how recently each cell spiked and gates these updates. In
three-cell PING, E-to-E coupling changes the E-assembly lag and
frequency; once STDP is active, that timing and the weights coevolve.

## Essential model

For $\Delta t=t_{\rm post}-t_{\rm pre}$, a representative pairwise rule is

$$
\Delta w = \begin{cases}
A_+\,e^{-\Delta t/\tau_+} & \Delta t>0\ \text{(post after pre)}\\[2pt]
-A_-\,e^{\Delta t/\tau_-} & \Delta t<0\ \text{(pre after post)}
\end{cases}
$$

The network implementation approximates this with a smooth,
voltage-triggered version: each E-cell carries an adaptation-like trace,
and every E-to-E weight is nudged up or down whenever its pre- or
postsynaptic cell crosses spike threshold, softly clamped between 0 and a
fixed upper bound `B`.

## Code examples

All eight examples live in one notebook, [`chapter40.ipynb`](chapter40.ipynb).

`simulate_abbott_song`/`plot_abbott_song` cover `ABBOTT_SONG`: the
piecewise-exponential timing rule $F_0$ and a windowed version $F$ that
tapers to zero as $|z|\to0$, with `tau_plus` and `K_plus` exposed through
sliders.

`simulate_rtm_voltage_trace_with_a`/`plot_rtm_voltage_trace_with_a` cover
`RTM_VOLTAGE_TRACE_WITH_A`: an RTM voltage trace together with an
adaptation-like variable $a$ that relaxes toward 0 while the cell is
depolarized, with `i_ext` and `C` exposed through sliders.

A shared block of RTM/WB gating functions, `tau_peak_function`/
`tau_d_q_function` (synapse rise-time solver), `derivative_three_cell_ping`,
`spike_detection`, and `simulate_three_cell_ping`/`reciprocal_g_ee` covers
the fixed-coupling two-E/one-I network reused by four examples:

- `THREE_CELL_PING_1` (baseline, no E-to-E coupling) and
  `THREE_CELL_PING_2` (weak vs. strong reciprocal coupling) are plotted
  with `plot_three_cell_raster`, and both have interactive sliders over
  the coupling strength.
- `sweep_three_cell_ping_ee`/`plot_ee_sweep` cover `THREE_CELL_PING_3`
  (reciprocal coupling sweep) and `THREE_CELL_PING_4` (one-way coupling
  sweep), each reporting the E1-E2 lag and E2 frequency as a function of
  coupling strength.

`simulate_three_cell_ping_5`/`plot_three_cell_ping_5` cover
`THREE_CELL_PING_5`: the same two-E/one-I network, but with the
reciprocal E-to-E weights evolving under STDP, tracked via `g_12`, `g_21`
and the shrinking E1-E2 `lags`. `g_ee0` is exposed through a slider. Its
50000-step explicit-Heun loop (`_three_cell_ping5_loop`) is
`@njit`-compiled.

`simulate_ping_with_stdp`/`plot_ping_with_stdp_raster`/
`plot_ping_with_stdp_density` cover `PING_WITH_STDP`: 200 RTM E-cells and
50 WB I-cells with every E-to-E synapse plastic under STDP, reporting the
raster and the final kernel-density estimate of E-to-E synaptic strengths
(`vec_g`). Its per-timestep update (`_stdp_loop`, shared with
`THREE_CELL_PING_5`'s STDP increment kernel `_g_ee_derivative_s`) is
`@njit`-compiled; a full run takes a few minutes even compiled, so no
interactive slider is provided for it.

## What to look for

Read the Abbott--Song curve with signed pre/post timing before
interpreting weight traces. In the fixed-coupling three-cell cases,
compare E-cell lag and frequency as coupling changes. In the STDP cases,
ask whether weight evolution changes the spike order that produced it:
in `THREE_CELL_PING_5`, the faster-driven E-cell reliably leads, so its
outgoing synapse potentiates toward its bound while the reverse synapse
decays, and the lag itself shrinks as that asymmetry grows. In
`PING_WITH_STDP`, look for a widened, non-degenerate spread in the final
E-to-E weight distribution rather than every synapse sitting at its
initial value.

## Suggested order

1. Run `ABBOTT_SONG` and `RTM_VOLTAGE_TRACE_WITH_A`.
2. Run `THREE_CELL_PING_1` through `THREE_CELL_PING_4`.
3. Finish with `THREE_CELL_PING_5` and `PING_WITH_STDP`.

## Prerequisites and related chapters

Chapter 9 provides adaptation context, Chapter 20 chemical synapses,
Chapter 30 PING, and Chapter 39 short-term plasticity.

## Running the examples

Open [`chapter40.ipynb`](chapter40.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom.
NumPy, SciPy, Matplotlib, ipywidgets, and numba are required. The
`THREE_CELL_PING_3`/`_4` coupling sweeps and the `PING_WITH_STDP` example
are the slowest cells and can take several minutes; that is expected.
