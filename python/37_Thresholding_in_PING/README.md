# Thresholding in PING

## Overview

PING can have sharp boundaries between suppression and participation. These examples isolate a non-reset threshold mechanism, show a baseline PING threshold calculation, and magnify the transition where a small timing or drive change changes cycle recruitment.

## Core ideas

Thresholding is a network property: excitation must arrive during the interval left by recurrent inhibition. A non-reset reference separates continuous voltage crossing from an imposed reset. Near a boundary, a cell can miss a cycle or recruit I cells into the next one.

## Essential model

The E-to-I and I-to-E PING conductance loop evolves continuously, while a voltage threshold selects spike events. In the non-reset comparison, crossing does not reset the voltage. The key output is the timing or parameter boundary between firing and nonparticipation.

## Code examples

- [`NO_RESET`](NO_RESET/) simulates the threshold mechanism without voltage reset.
- [`PING_THR_1`](PING_THR_1/) shows a baseline PING thresholding calculation.
- [`PING_THR_1_ZOOM`](PING_THR_1_ZOOM/) magnifies the suppression/participation transition.
- [`THRESHOLDING`](THRESHOLDING/) plots the thresholding construction directly.

## What to look for

Use `NO_RESET` to distinguish a threshold event from reset dynamics. In the zoom, inspect whether E crosses when inhibition relaxes enough to recruit I, or remains suppressed. A boundary represents changed cycle participation, not merely a small voltage difference.

## Suggested order

1. Run `NO_RESET` and `THRESHOLDING`.
2. Run `PING_THR_1`.
3. Use `PING_THR_1_ZOOM` to inspect its boundary.

## Prerequisites and related chapters

Chapter 7 discusses reset dynamics, Chapter 30 PING, Chapter 35 inhibitory windows, and Chapter 38 gamma coherence.

## Running the examples

Run `python main.py` from the selected immediate example directory. NumPy and Matplotlib are required.
