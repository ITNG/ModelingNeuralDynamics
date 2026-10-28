# ModelingNeuralDynamics
An Introduction to Modeling Neuronal Dynamics - Christoph Borgers in python

[![tests](https://github.com/ITNG/ModelingNeuralDynamics/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/ITNG/ModelingNeuralDynamics/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/modelingneuraldynamics.svg)](https://pypi.org/project/modelingneuraldynamics/)
[![Python versions](https://img.shields.io/pypi/pyversions/modelingneuraldynamics.svg)](https://pypi.org/project/modelingneuraldynamics/)
[![License: GPL v3](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)

<p align="center">
<img src="https://raw.githubusercontent.com/ITNG/ModelingNeuralDynamics/main/docs/img/chapter30_ping4.png" width="32%">
<img src="https://raw.githubusercontent.com/ITNG/ModelingNeuralDynamics/main/docs/img/chapter38_gamma_coherence_1.png" width="32%">
<img src="https://raw.githubusercontent.com/ITNG/ModelingNeuralDynamics/main/docs/img/chapter40_ping_with_stdp.png" width="32%">
</p>
<p align="center">
PING gamma rhythm &nbsp;&middot;&nbsp; gamma coherence &nbsp;&middot;&nbsp; PING network with STDP
</p>

### Installation

The shared helper package used by some chapters is on PyPI:

```bash
pip install modelingneuraldynamics
```

For a local checkout, create a uv environment and install this package with
the notebook and development dependencies:

```bash
uv sync --all-groups
uv run jupyter lab
```

Run Brian2 regression tests explicitly with
`uv run pytest -m 'brian2 and not slow'`; the default test selection excludes them.

After cloning, install the notebook Git filter:

```bash
uv run python scripts/strip_notebook.py --install
```

The filter keeps notebook and cell metadata empty, clears outputs and execution
counts, and normalizes cell IDs when staging. Local notebooks retain their
metadata and outputs. Git configuration is local, so run setup in each clone.

### What's supported

Chapters are ported from the book's original MATLAB programs to Python.
Every tracked chapter is now a single, tested `chapterNN.ipynb` notebook
that installs its own dependencies and opens directly in Colab.

- **38 / 38** tracked chapters converted to notebooks (✅)
- Chapters 2 and 6 have no Python example in the book.

<details>
<summary>Full chapter-by-chapter status (click to expand)</summary>

| # | Chapter | Status | Open |
|---|---|---|---|
| 1 | Modeling a Single Neuron | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter01.ipynb) [guide](python/chapter01.md) |
| 3 | The Classical HH ODEs | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter03.ipynb) [guide](python/chapter03.md) |
| 4 | Numerical Solution of HH ODEs | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter04.ipynb) [guide](python/chapter04.md) |
| 5 | The Simple Model of Neurons in Rodent Brains | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter05.ipynb) [guide](python/chapter05.md) |
| 7 | Linear Integrate-and-Fire (LIF) Neurons | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter07.ipynb) [guide](python/chapter07.md) |
| 8 | Quadratic Integrate-and-Fire (QIF) and Theta Neurons | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter08.ipynb) [guide](python/chapter08.md) |
| 9 | Spike Frequency Adaptation | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter09.ipynb) [guide](python/chapter09.md) |
| 10 | The Slow-Fast Phase Plane | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter10.ipynb) [guide](python/chapter10.md) |
| 11 | The Saddle-Node Bifurcation | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter11.ipynb) [guide](python/chapter11.md) |
| 12 | Two-Dimensional Bifurcation Analysis | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter12.ipynb) [guide](python/chapter12.md) |
| 13 | Hopf Bifurcations | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter13.ipynb) [guide](python/chapter13.md) |
| 14 | Model Neurons of Bifurcation Type 2 | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter14.ipynb) [guide](python/chapter14.md) |
| 15 | Canard Explosions | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter15.ipynb) [guide](python/chapter15.md) |
| 16 | Model Neurons of Bifurcation Type 3 | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter16.ipynb) [guide](python/chapter16.md) |
| 17 | Frequency-Current Curves | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter17.ipynb) [guide](python/chapter17.md) |
| 18 | Bistability Resulting from Rebound Firing | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter18.ipynb) [guide](python/chapter18.md) |
| 19 | Bursting | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter19.ipynb) [guide](python/chapter19.md) |
| 20 | Chemical Synapses | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter20.ipynb) [guide](python/chapter20.md) |
| 21 | Gap Junctions | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter21.ipynb) [guide](python/chapter21.md) |
| 22 | A Wilson-Cowan Model of an Oscillatory E-I Network | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter22.ipynb) [guide](python/chapter22.md) |
| 23 | Entrainment by Excitatory Input Pulses | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter23.ipynb) [guide](python/chapter23.md) |
| 24 | Synchronization by Fast Recurrent Excitation | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter24.ipynb) [guide](python/chapter24.md) |
| 25 | Phase Response Curves (PRCs) | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter25.ipynb) [guide](python/chapter25.md) |
| 26 | Phase Locking of Two Oscillators | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter26.ipynb) [guide](python/chapter26.md) |
| 27 | Phase Locking with Delays | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter27.ipynb) [guide](python/chapter27.md) |
| 28 | Weakly Coupled Oscillators | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter28.ipynb) [guide](python/chapter28.md) |
| 29 | Stability of the Synchronous State | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter29.ipynb) [guide](python/chapter29.md) |
| 30 | The PING Model of Gamma Rhythms | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter30.ipynb) [guide](python/chapter30.md) |
| 31 | ING Rhythms | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter31.ipynb) [guide](python/chapter31.md) |
| 32 | M-Current PING and Poisson PING | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter32.ipynb) [guide](python/chapter32.md) |
| 33 | M-Current PING and PINB | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter33.ipynb) [guide](python/chapter33.md) |
| 34 | Nested Gamma-Theta Rhythms | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter34.ipynb) [guide](python/chapter34.md) |
| 35 | Periodic Inhibition | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter35.ipynb) [guide](python/chapter35.md) |
| 36 | F-I Curves: Pulsed Excitation | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter36.ipynb) [guide](python/chapter36.md) |
| 37 | Thresholding in PING | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter37.ipynb) [guide](python/chapter37.md) |
| 38 | Gamma Coherence | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter38.ipynb) [guide](python/chapter38.md) |
| 39 | Short-Term Depression and Facilitation | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter39.ipynb) [guide](python/chapter39.md) |
| 40 | Spike-Timing-Dependent Plasticity (STDP) | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter40.ipynb) [guide](python/chapter40.md) |

</details>

### Running the Brian2 chapters

`brian/` contains **42 notebooks covering 37 chapters**, including representative
examples for 23 newly supported chapters. Shared neuron, phase, network, forcing,
and plasticity modules live in `mnd/brian/`.

See the [Brian coverage and setup guide](brian/README.md) for every notebook,
its implemented examples, and the distinction between representative coverage
and complete Python figure parity. Chapter 03 is retained in Python because its
examples are entirely static gating-curve calculations; chapters 02 and 06 have
no Python examples.

New Brian notebooks require this source checkout and its uv environment; the new
helpers are not assumed to be present in an older PyPI release.

### Introduction 
This book is intended as a text for a one-semester course on Mathematical and Computational Neuroscience for upper-level undergraduate and beginning graduate students of mathematics, the natural sciences, engineering, or computer science. An undergraduate introduction to differential equations is more than enough mathematical background. Only a slim, high school-level background in physics is assumed, and none in biology.

Topics include models of individual nerve cells and their dynamics, models of networks of neurons coupled by synapses and gap junctions, origins and functions of population rhythms in neuronal networks, and models of synaptic plasticity.

An extensive online collection of Matlab programs generating the figures accompanies the book.

### matlab code gathered from [here](https://link.springer.com/book/10.1007/978-3-319-51171-9)


### Python codes provided by contributors
See the [practical Python chapter guides](python/README.md) for the concepts,
equations, example map, and expected results for every implemented chapter.
