# Chapter 35: Periodic Inhibition

## Overview

Periodic inhibitory forcing creates recurring response windows rather than merely lowering firing. The examples show the shape of the periodic gate, oscillatory LIF responses (deterministic and noisy) under periodic vs. tonic inhibition, and f-I curves under periodic inhibition for both the LIF and the RTM neuron.

## Core ideas

Spikes can be locked, suppressed, or shifted according to their phase in the inhibitory cycle. Inhibited f-I curves count spikes that survive these windows, so changes in average rate can reflect cycle skipping rather than a uniformly slower intrinsic response.

## Essential model

The forcing is $I_{\rm inh}(t)=g_{\rm inh}s(t)(E_I-v)$, with a periodic gate $s(t) \propto \exp(\alpha\cos^2(\pi t/T)) - 1$, normalized so its average over one period equals a chosen mean conductance $\bar g$. Sweeping applied current and counting spikes per observation time under this forcing produces the inhibited f-I relation.

## Code examples

All examples live in [`chapter35.ipynb`](chapter35.ipynb).

- `phi_of` / `plot_oscillations` / `plot_phi` (`OSCILLATIONS`): the normalized periodic gate shape for different sharpness values $\alpha$.
- `run_periodic_inhibition` / `plot_periodic_inhibition` (`PERIODIC_INHIBITION`, `PERIODIC_INHIBITION_3`): LIF neuron under periodic vs. tonic (mean) inhibitory conductance, for gate sharpness $\alpha=5$ and $\alpha=1$ respectively.
- `run_periodic_inhibition_noisy` / `plot_periodic_inhibition_noisy` (`PERIODIC_INHIBITION_2`): the same comparison with an added Ornstein-Uhlenbeck noise current.
- `compute_periodic_inhibition_f_i_curve` / `compute_periodic_inhibition_f_i_curve_2` / `plot_periodic_inhibition_f_i_curve` (`PERIODIC_INHIBITION_F_I_CURVE`, `PERIODIC_INHIBITION_F_I_CURVE_2`): step-like f-I curves of the LIF neuron under periodic inhibition ($\alpha=5$ and $\alpha=1$), compared against the closed-form tonic-inhibition f-I curve. The sweep's inner loop is `@njit`-compiled.
- `compute_rtm_f_i_curves` / `compute_rtm_f_i_curves_2` / `plot_rtm_f_i_curves` (`RTM_F_I_CURVE_WITH_INHIBITION`, `RTM_F_I_CURVE_WITH_INHIBITION_2`): tonic vs. periodic-inhibition f-I curves for the conductance-based RTM neuron, for periodic peak amplitude $\bar g$ and $2\bar g$. Both sweeps are `@njit`-compiled.

## What to look for

Follow spikes relative to inhibitory cycles, especially at a response-window boundary. Compare f-I slopes, offsets, and plateaus: lower rate can arise from skipped cycles.

## Suggested order

1. Run `plot_oscillations` (and the `plot_phi` widget) to see the gate shape.
2. Run `plot_periodic_inhibition` for $\alpha=5$ and $\alpha=1$, then `plot_periodic_inhibition_noisy`.
3. Compare the two periodic-inhibition LIF f-I curves.
4. Contrast the two RTM f-I calculations.

## Prerequisites and related chapters

Chapter 17 introduces steady f-I curves, Chapter 23 periodic forcing, Chapter 31 ING, and Chapter 36 pulsed excitation.

## Running the examples

Open [`chapter35.ipynb`](chapter35.ipynb) and run the cells in order (or use Colab via the badge at the top of the notebook). NumPy, Matplotlib, ipywidgets, and Numba are required; the f-I sweeps take longer than a single trace.
