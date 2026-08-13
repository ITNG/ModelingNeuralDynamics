# Python chapter guides

This directory contains Python programs that reproduce and extend the
computational examples in *An Introduction to Modeling Neuronal Dynamics*.
The guides below explain the concepts, equations, examples, and expected
results behind each implemented chapter.

Repository numbering follows the historical code organization, so it is not
always an exact mirror of the book's table of contents. Chapters 2 and 6 have
no Python example directories. The Chapter 9 guide is deferred and is not
linked from this index.

## Foundations and single-neuron models (Chapters 1-9)

- [Chapter 1: Modeling a Single Neuron](chapter01.md)
- [Chapter 3: The Classical HH ODEs](chapter03.md)
- [Chapter 4: Numerical Solution of HH ODEs](chapter04.md)
- [Chapter 5: The Simple Model of Neurons in Rodent Brains](chapter05.md)
- [Chapter 7: Linear Integrate-and-Fire (LIF) Neurons](chapter07.md)
- [Chapter 8: Quadratic Integrate-and-Fire (QIF) and Theta Neurons](chapter08.md)

## Single-neuron dynamics (Chapters 10-19)

- [Chapter 10: The Slow-Fast Phase Plane](chapter10.md)
- [Chapter 11: The Saddle-Node Bifurcation](chapter11.md)
- [Chapter 12: Two-Dimensional Bifurcation Analysis](chapter12.md)
- [Chapter 13: Hopf Bifurcations](chapter13.md)
- [Chapter 14: Model Neurons of Bifurcation Type 2](chapter14.md)
- [Chapter 15: Canard Explosions](chapter15.md)
- [Chapter 16: Model Neurons of Bifurcation Type 3](chapter16.md)
- [Chapter 17: Frequency-Current Curves](chapter17.md)
- [Chapter 18: Bistability Resulting from Rebound Firing](chapter18.md)
- [Chapter 19: Bursting](chapter19.md)

## Communication (Chapters 20-22)

- [Chapter 20: Chemical Synapses](chapter20.md)
- [Chapter 21: Gap Junctions](chapter21.md)
- [Chapter 22: A Wilson-Cowan Model of an Oscillatory E-I Network](chapter22.md)

## Entrainment, synchronization, and rhythms (Chapters 23-38)

- [Chapter 23: Entrainment by Excitatory Input Pulses](chapter23.md)
- [Chapter 24: Synchronization by Fast Recurrent Excitation](chapter24.md)
- [Chapter 25: Phase Response Curves (PRCs)](chapter25.md)
- [Chapter 26: Phase Locking of Two Oscillators](chapter26.md)
- [Chapter 27: Phase Locking with Delays](chapter27.md)
- [Chapter 28: Weakly Coupled Oscillators](chapter28.md)
- [Chapter 29: Stability of the Synchronous State](chapter29.md)
- [Chapter 30: The PING Model of Gamma Rhythms](chapter30.md)
- [Chapter 31: ING Rhythms](chapter31.md)
- [Chapter 32: M-Current PING and Poisson PING](chapter32.md)
- [Chapter 33: M-Current PING and PINB](chapter33.md)
- [Chapter 34: Nested Gamma-Theta Rhythms](chapter34.md)
- [Chapter 35: Periodic Inhibition](chapter35.md)
- [Chapter 36: F-I Curves: Pulsed Excitation](chapter36.md)
- [Chapter 37: Thresholding in PING](chapter37.md)
- [Chapter 38: Gamma Coherence](chapter38.md)

## Plasticity (Chapters 39-40)

- [Chapter 39: Short-Term Depression and Facilitation](chapter39.md)
- [Chapter 40: Spike-Timing-Dependent Plasticity (STDP)](chapter40.md)

## How to use these guides

Start with each guide's overview and core ideas, then use the essential model
to connect equations to their symbols. The code-example map identifies the
programs that implement each result, while the interpretation section explains
what to look for in the output. Finish with the run instructions to reproduce
or adapt an example locally.
