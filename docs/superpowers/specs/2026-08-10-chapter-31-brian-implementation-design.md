# Chapter 31 Brian2 Implementation Design

## Goal

Evaluate every Chapter 31 Python example for genuine Brian2 compatibility and
implement the compatible neuronal simulations. Chapter 31 contains 16 Python
examples. Fourteen use simulated neuron dynamics and belong in Brian2; the two
abstract pulse-coupling calculations are analytical phase-map examples and
will be marked `n/a`.

The compatibility classification is:

- Brian2: `1_CELL_ING`, `1_CELL_ING_CONDITION_NUMBERS`, `ING_1` through
  `ING_10`, `ING_ENTRAINING_E_CELLS`, and `ING_ENTRAINING_E_CELLS_2`.
- Not applicable: `ABSTRACT_PULSE_COUPLING_INH` and
  `ABSTRACT_PULSE_COUPLING_INH_2`.

The MATLAB implementations remain the source reference. The existing Python
ports, which already have MATLAB-backed tests, are the primary numeric oracle
for the Brian2 implementation.

## Structure

Create `brian/chapter31.ipynb` as the single Chapter 31 Brian2 notebook. Keep
simulation functions independent of plotting cells so tests can load notebook
definitions without running the full examples.

Organize the notebook into four layers:

1. Shared Wang-Buzsaki inhibitory-neuron equations, synaptic `q`/`s`
   dynamics, initialization helpers, deterministic connectivity generation,
   and raster plotting.
2. `simulate_single_cell_ing(...)` and
   `compute_ing_condition_numbers(...)` for the two single-cell examples.
3. `simulate_ing_network(config, ...)`, driven by an `ING_CONFIGS` mapping for
   `ING_1` through `ING_10`.
4. `simulate_ing_entrainment(...)` for the two networks combining RTM
   excitatory neurons and WB inhibitory neurons.

Simulation functions return explicit result objects containing spike times,
neuron indices, relevant state traces, and the effective configuration.
Plotting cells consume these results to reproduce the reference traces,
rastergrams, condition-number table, and entrainment comparisons.

## Brian2 Models

The inhibitory population uses the WB membrane and gating equations from the
Python reference, with instantaneous sodium activation and dynamic potassium
activation and sodium inactivation. Continuous synaptic release follows the
reference `q` and `s` equations.

Chemical inhibition and electrical coupling use Brian2 `Synapses`. Chemical
current is summed from presynaptic `s` values. Gap-junction current is a
bidirectional diffusive current proportional to the voltage difference. The
implementation preserves each example's self-connection convention,
fixed-versus-random indegree, normalization, and coupling probabilities.

The entrainment model adds an RTM excitatory population with its own membrane
and gating equations and uses separate E-to-E, E-to-I, I-to-E, I-to-I, and
gap-junction connection sets matching the Python references.

The condition-number example remains Brian2-compatible because its derived
sensitivities come from repeated neuronal simulations. Only the final
percentage calculation is ordinary numerical postprocessing.

## Configuration and Data Flow

Each configuration records the corresponding Python values for population
sizes, heterogeneous external drive, chemical and electrical connection
strengths and probabilities, synaptic time constants, fixed-versus-random
indegree, duration, timestep, and seed.

The execution flow is:

1. A NumPy random generator creates external currents, adjacency matrices, and
   initial-state samples deterministically.
2. The generated arrays are assigned directly to Brian2 neuron groups and
   synapses. This avoids accidental differences between NumPy and Brian random
   streams.
3. Brian2 integrates the model with the reference timestep and an appropriate
   explicit integration method.
4. Spike and state monitors are converted into stable NumPy-facing results.
5. Plotting cells render the example-specific output from those results.

Every public simulation starts with `brian2.start_scope()` so repeated tests
and notebook executions do not retain objects from earlier runs.

## Validation and Error Handling

Public helpers validate population sizes, probabilities, coupling strengths,
and optional array shapes. Probabilities must lie in `[0, 1]`. A nonzero
coupling strength with zero connection probability is rejected because its
normalized weight is undefined. A zero coupling strength creates no synapses
and never performs a division by zero.

Validation follows the existing provenance chain:

`MATLAB reference -> verified Python port -> Brian2 port`

Tests are added family-by-family:

- Single-cell ING tests compare the voltage behavior, spike period, and firing
  frequency with the Python port.
- The condition-number test compares the three percentage sensitivities for
  input current, inhibitory conductance, and synaptic decay.
- Tests for `ING_1` through `ING_10` verify every configuration and compare
  activity level, active-cell fraction, population rate, and binned raster
  structure. Identical generated currents and connectivity are supplied to
  Python and Brian2 where practical.
- Entrainment tests compare E/I spike counts, active-cell fractions,
  population frequencies, and response ordering across the three drive
  conditions.
- A notebook smoke test loads all definitions without executing plotting
  cells.

Focused tests use shortened deterministic runs when those runs still detect
incorrect equations, normalization, connectivity, or parameter mappings.
Final verification also runs every example at its full reference duration and
regenerates figures for visual comparison with the Python output.

## Completion Criteria

Chapter 31 is complete when:

- all 14 compatible examples have callable Brian2 implementations and their
  corresponding notebook presentation cells;
- all focused and Chapter 31 aggregate tests pass;
- full-duration simulations complete and their figures have been visually
  checked against the Python references;
- `PROGRESS.md` marks the 14 simulations complete and the two analytical
  pulse-coupling examples `n/a`; and
- no verified Python or MATLAB implementation has been changed.

Existing unrelated work in Chapters 19, 20, 34, and 35 must not be modified or
staged as part of this work.
