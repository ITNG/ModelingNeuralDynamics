# The saddle-node bifurcation

## Overview

This chapter's single example makes a saddle-node collision concrete: in the
displayed positive-$x$ region, a stable fixed point and a saddle coalesce
into one half-stable point and then disappear.

## Core ideas

A saddle-node bifurcation changes the number of equilibria in a local region.
Before the displayed collision, trajectories can be attracted to a node or
separated by a saddle's stable manifold. At the collision, the flow is
attracting from one side and repelling from the other; beyond it, the displayed
positive-$x$ pair no longer remains.

## Essential model

The plotted planar normal-form-like system is

$$
\dot x=-ax+y,\qquad \dot y=\frac{x^2}{1+x^2}-by.
$$

Here $x$ and $y$ are state variables, $a$ is the varied control
parameter, and $b=1$ is fixed. In the displayed positive-$x$ region,
solving the nullclines produces two fixed points for $a<0.5$, one collision
point at $a=0.5$, and none for $a>0.5$. The full system also has the
equilibrium $(x,y)=(0,0)$ for every $a$; the script clips it out by
starting and displaying trajectories at $x\ge0.2$.

## Code examples

The example now lives in [`chapter11.ipynb`](chapter11.ipynb):
`plot_saddle_node_bifurcation` integrates trajectories for the three
parameter regimes and draws their phase portraits, using
`saddle_node_fixed_points` for the closed-form fixed points.

## What to look for

Read the three panels from left to right. The black filled marker is the
stable fixed point and the white marker is the saddle; their merged half-filled
marker in the middle panel shows the collision. In the final panel, all shown
trajectories pass through because no positive-$x$ fixed point survives; the
equilibrium at the excluded origin is outside the plotted frame.

## Suggested order

1. Run `SADDLE_NODE_BIFURCATION` and identify the fixed-point markers.
2. Compare nearby trajectories on each side of the saddle's separatrix.
3. Relate the three panels to creation or destruction of a resting state.

## Prerequisites and related chapters

Chapter 10 supplies nullclines and phase portraits. Chapter 12 applies
saddle-node ideas to reduced neuron dynamics, while Chapters 17--18 connect
them to firing thresholds and bistability.

## Running the examples

Open [`chapter11.ipynb`](chapter11.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook, and run all cells top to bottom.
