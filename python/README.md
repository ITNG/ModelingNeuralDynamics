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
- [Chapter 25: Phase Response Curves (PRCs)](25_Phase_Response_Curves_%28PRCs%29/README.md)
- [Chapter 26: Phase Locking of Two Oscillators](26_Phase_Locking_of_Two_Oscillators/README.md)
- [Chapter 27: Phase Locking with Delays](27_Phase_Locking_with_Delays/README.md)
- [Chapter 28: Weakly Coupled Oscillators](28_Weakly_Coupled_Oscillators/README.md)
- [Chapter 29: Stability of the Synchronous State](29_Stability_of_the_Synchronous_State/README.md)
- [Chapter 30: The PING Model of Gamma Rhythms](30_The_PING_Model_of_Gamma_Rhythms/README.md)
- [Chapter 31: ING Rhythms](31_ING_Rhythms/README.md)
- [Chapter 32: M-Current PING and Poisson PING](32_M_Current_PING_and_Poisson_PING/README.md)
- [Chapter 33: M-Current PING and PINB](33_M_Current_PING_and_PINB/README.md)
- [Chapter 34: Nested Gamma-Theta Rhythms](34_Nested_Gamma_Theta_Rhythms/README.md)
- [Chapter 35: Periodic Inhibition](35_Periodic_Inhibition/README.md)
- [Chapter 36: F-I Curves: Pulsed Excitation](36_F_I_Curves_Pulsed_Excitation/README.md)
- [Chapter 37: Thresholding in PING](37_Thresholding_in_PING/README.md)
- [Chapter 38: Gamma Coherence](38_Gamma_Coherence/README.md)

## Plasticity (Chapters 39-40)

- [Chapter 39: Short-Term Depression and Facilitation](39_Short-Term_Depression_and_Facilitation/README.md)
- [Chapter 40: Spike-Timing-Dependent Plasticity (STDP)](40_Spike_Timing-Dependent_Plasticity%28STDP%29/README.md)

## How to use these guides

Start with each guide's overview and core ideas, then use the essential model
to connect equations to their symbols. The code-example map identifies the
programs that implement each result, while the interpretation section explains
what to look for in the output. Finish with the run instructions to reproduce
or adapt an example locally.
