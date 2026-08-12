# ModelingNeuralDynamics
An Introduction to Modeling Neuronal Dynamics - Christoph Borgers in python

[![tests](https://github.com/ITNG/ModelingNeuralDynamics/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/ITNG/ModelingNeuralDynamics/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/modelingneuraldynamics.svg)](https://pypi.org/project/modelingneuraldynamics/)
[![Python versions](https://img.shields.io/pypi/pyversions/modelingneuraldynamics.svg)](https://pypi.org/project/modelingneuraldynamics/)
[![License: GPL v3](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)

<p align="center">
<img src="https://raw.githubusercontent.com/ITNG/ModelingNeuralDynamics/main/python/30_The_PING_Model_of_Gamma_Rhythms/PING_4/fig.png" width="32%">
<img src="https://raw.githubusercontent.com/ITNG/ModelingNeuralDynamics/main/python/38_Gamma_Coherence/GAMMA_COHERENCE_1/fig.png" width="32%">
<img src="https://raw.githubusercontent.com/ITNG/ModelingNeuralDynamics/main/python/40_Spike_Timing-Dependent_Plasticity%28STDP%29/PING_WITH_STDP/fig1.png" width="32%">
</p>
<p align="center">
PING gamma rhythm &nbsp;&middot;&nbsp; gamma coherence &nbsp;&middot;&nbsp; PING network with STDP
</p>

### Installation

The shared helper package used by some chapters is on PyPI:

```bash
pip install modelingneuraldynamics
```

### What's supported

Chapters are being ported from the book's original MATLAB programs to
Python one at a time. A converted chapter is a single, tested
`chapterNN.ipynb` notebook that installs its own dependencies and opens
directly in Colab; chapters not yet converted still live as one
`main.py` script per example under `python/<NN_chapter_name>/`.

- **22 / 37** tracked chapters converted to notebooks (✅)
- **15** chapters still as legacy per-example scripts (🚧)
- Chapters 2 and 6 have no Python example in the book. Chapter 9 has a
  pilot notebook (below) but isn't part of the official guide index yet.

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
| 9 | Spike Frequency Adaptation | ✅ Notebook | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/python/chapter09.ipynb) guide pending |
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
| 26 | Phase Locking of Two Oscillators | 🚧 Legacy scripts | [guide](python/26_Phase_Locking_of_Two_Oscillators/README.md) |
| 27 | Phase Locking with Delays | 🚧 Legacy scripts | [guide](python/27_Phase_Locking_with_Delays/README.md) |
| 28 | Weakly Coupled Oscillators | 🚧 Legacy scripts | [guide](python/28_Weakly_Coupled_Oscillators/README.md) |
| 29 | Stability of the Synchronous State | 🚧 Legacy scripts | [guide](python/29_Stability_of_the_Synchronous_State/README.md) |
| 30 | The PING Model of Gamma Rhythms | 🚧 Legacy scripts | [guide](python/30_The_PING_Model_of_Gamma_Rhythms/README.md) |
| 31 | ING Rhythms | 🚧 Legacy scripts | [guide](python/31_ING_Rhythms/README.md) |
| 32 | M-Current PING and Poisson PING | 🚧 Legacy scripts | [guide](python/32_M_Current_PING_and_Poisson_PING/README.md) |
| 33 | M-Current PING and PINB | 🚧 Legacy scripts | [guide](python/33_M_Current_PING_and_PINB/README.md) |
| 34 | Nested Gamma-Theta Rhythms | 🚧 Legacy scripts | [guide](python/34_Nested_Gamma_Theta_Rhythms/README.md) |
| 35 | Periodic Inhibition | 🚧 Legacy scripts | [guide](python/35_Periodic_Inhibition/README.md) |
| 36 | F-I Curves: Pulsed Excitation | 🚧 Legacy scripts | [guide](python/36_F_I_Curves_Pulsed_Excitation/README.md) |
| 37 | Thresholding in PING | 🚧 Legacy scripts | [guide](python/37_Thresholding_in_PING/README.md) |
| 38 | Gamma Coherence | 🚧 Legacy scripts | [guide](python/38_Gamma_Coherence/README.md) |
| 39 | Short-Term Depression and Facilitation | 🚧 Legacy scripts | [guide](python/39_Short-Term_Depression_and_Facilitation/README.md) |
| 40 | Spike-Timing-Dependent Plasticity (STDP) | 🚧 Legacy scripts | [guide](python/40_Spike_Timing-Dependent_Plasticity%28STDP%29/README.md) |

</details>

### Running `brian/` chapters on Colab

`brian/` holds a separate, Brian2-based implementation, one notebook per
chapter. Every notebook there can also be opened directly in Google
Colab — click a chapter's badge below.

| Chapter | Colab |
|---|---|
| 01 - Modeling a Single Neuron | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter01.ipynb) |
| 04 - Numerical Solution of HH ODEs | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter04.ipynb) |
| 05 - Three Simple Models of Neurons in Rodent Brains | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter05.ipynb) |
| 07 - Linear Integrate and Fire (LIF) Neurons | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter07.ipynb) |
| 08 - Quadratic Integrate and Fire (QIF) and Theta Neurons | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter08.ipynb) |
| 09 - Spike Frequency Adaptation | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter09.ipynb) |
| 20 - Chemical Synapses | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter20.ipynb) |

### Introduction 
This book is intended as a text for a one-semester course on Mathematical and Computational Neuroscience for upper-level undergraduate and beginning graduate students of mathematics, the natural sciences, engineering, or computer science. An undergraduate introduction to differential equations is more than enough mathematical background. Only a slim, high school-level background in physics is assumed, and none in biology.

Topics include models of individual nerve cells and their dynamics, models of networks of neurons coupled by synapses and gap junctions, origins and functions of population rhythms in neuronal networks, and models of synaptic plasticity.

An extensive online collection of Matlab programs generating the figures accompanies the book.

### matlab code gathered from [here](https://link.springer.com/book/10.1007/978-3-319-51171-9)


### Python codes provided by contributors
See the [practical Python chapter guides](python/README.md) for the concepts,
equations, example map, and expected results for every implemented chapter.
