# ING rhythms

## Overview

Interneuron gamma (ING) is an inhibitory network rhythm in which recovery from shared inhibition determines the next spike window. These examples include one-cell references, inhibitory pulse maps, population clustering and optional electrical coupling, plus ING entrainment of E cells.

## Core ideas

ING synchrony is organized by inhibitory recovery rather than recurrent excitation. Inhibitory phase-reset maps predict stable timing; gap junctions can alter clustering. The population outputs show whether an inhibitory rhythm creates repeatable windows for excitatory spikes.

## Essential model

Recurrent inhibition has the form $I_{{\rm II},i}=g_{\rm II}\sum_j s_j(E_I-v_i)$, optionally supplemented by gap-junction voltage-difference currents. A pulse map $\phi\mapsto F(\phi)$ reduces inhibitory timing to fixed points and their slopes.

## Code examples

- [`1_CELL_ING`](1_CELL_ING/) gives a single inhibitory-cell timing reference.
- [`1_CELL_ING_CONDITION_NUMBERS`](1_CELL_ING_CONDITION_NUMBERS/) measures timing sensitivity.
- [`ABSTRACT_PULSE_COUPLING_INH`](ABSTRACT_PULSE_COUPLING_INH/) and [`ABSTRACT_PULSE_COUPLING_INH_2`](ABSTRACT_PULSE_COUPLING_INH_2/) plot inhibitory pulse-coupling maps.
- [`ING_1`](ING_1/) simulates a population and writes raster/LFP output.
- [`ING_2`](ING_2/), [`ING_3`](ING_3/), [`ING_4`](ING_4/), [`ING_5`](ING_5/), [`ING_6`](ING_6/), [`ING_7`](ING_7/), [`ING_8`](ING_8/), [`ING_9`](ING_9/), and [`ING_10`](ING_10/) compare inhibitory population configurations.
- [`ING_ENTRAINING_E_CELLS`](ING_ENTRAINING_E_CELLS/) and [`ING_ENTRAINING_E_CELLS_2`](ING_ENTRAINING_E_CELLS_2/) show E-cell entrainment by ING.

## What to look for

Cells should fire after inhibitory conductance decays enough to open a recovery window. Compare abstract-map fixed points with persistent raster timing. In the entrainment cases, identify E spikes relative to the inhibitory window rather than expecting every E cell to fire every cycle.

## Suggested order

1. Run the one-cell and abstract-map examples.
2. Compare `ING_1` through `ING_10`.
3. Finish with the two `ING_ENTRAINING_E_CELLS` cases.

## Prerequisites and related chapters

Chapter 21 introduces gap junctions, Chapters 25--29 phase maps, Chapter 30 PING, and Chapter 35 periodic inhibition.

## Running the examples

Run `python main.py` in an immediate example directory. `ING_1` writes spike-time and LFP files. Network runs use NumPy, SciPy, and Matplotlib and can be slower than the maps.
