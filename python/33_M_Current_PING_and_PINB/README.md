# Beta rhythms from M-current PING and PINB

## Overview

The repository's historical title, “M Current PING and PINB,” covers the book's beta-rhythm material. The examples examine M-current PING period skipping and PING-with-inhibition/beta (PINB), including gap junctions and cell-assembly timing.

## Core ideas

Slow M-current recovery can make E cells skip gamma opportunities, yielding beta-scale population timing. PINB organizes which E assemblies participate and when; gap junctions can change inhibitory coordination. A slower rhythm therefore need not mean that every cell fires slowly.

## Essential model

The PING E--I loop is supplemented by M-current-controlled E-cell availability. Population period and skipped cycles are useful timing measures. Gap-junction currents are proportional to $v_j-v_i$, providing electrical coupling without chemical-synapse delay.

## Code examples

- [`M_CURRENT_BETA_WITH_GJ`](M_CURRENT_BETA_WITH_GJ/) simulates M-current beta timing with gap junctions.
- [`M_CURRENT_PING_4`](M_CURRENT_PING_4/), [`M_CURRENT_PING_5`](M_CURRENT_PING_5/), [`M_CURRENT_PING_6`](M_CURRENT_PING_6/), [`M_CURRENT_PING_7`](M_CURRENT_PING_7/), and [`M_CURRENT_PING_8`](M_CURRENT_PING_8/) compare period-skipping M-current PING cases.
- [`PINB_1`](PINB_1/), [`PINB_2`](PINB_2/), and [`PINB_3`](PINB_3/) simulate PINB cell-assembly configurations.

## What to look for

Count skipped gamma opportunities before calling the population rhythm beta. In the gap-junction case, inspect cluster timing, not only voltage. In PINB rasters, identify participating E assemblies and their order relative to inhibition.

## Suggested order

1. Run `M_CURRENT_PING_4` through `M_CURRENT_PING_8`.
2. Compare `M_CURRENT_BETA_WITH_GJ`.
3. Run `PINB_1`, `PINB_2`, and `PINB_3`.

## Prerequisites and related chapters

Chapter 9 introduces M-current context, Chapter 21 gap junctions, Chapters 30 and 32 PING, and Chapter 34 nested slow/fast rhythms.

## Running the examples

Run `python main.py` in an immediate example directory. The network simulations require NumPy, SciPy, and Matplotlib.
