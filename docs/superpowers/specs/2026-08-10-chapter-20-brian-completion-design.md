# Chapter 20 Brian2 Completion Design

## Goal

Complete the four remaining Brian2-compatible Chapter 20 examples using the
already verified Python ports as the primary numeric references:

- `B_JAHR_STEVENS`
- `RTM_PLOT_Q`
- `RTM_PLOT_S_PRESCRIBE_TAU_PEAK`
- `RTM_WITH_AUTAPSE_F_I_CURVE`

Work remains on `develop_python`. Existing uncommitted work in chapters 19,
33, and 34 must not be modified or staged.

## Structure

Extend `brian/chapter20.ipynb`, which already contains the other four Chapter
20 Brian2 examples. Add one clearly labeled section for each missing example.
Keep simulation and calculation logic in small callable functions so tests can
load notebook definitions without executing plotting cells.

Reuse the notebook's existing Reduced Traub-Miles equations where their
behavior and parameters match. Add narrowly scoped helpers only where an
example requires distinct behavior, such as prescribed synaptic peak timing or
autaptic current. Do not introduce a chapter-level package or unrelated shared
abstraction.

## Example Behavior

### B_JAHR_STEVENS

Provide a function that evaluates the voltage-dependent magnesium block over
the same voltage grid and with the same constants as the Python reference.
Add a plotting cell reproducing the reference curve. Although this example is
algebraic rather than a time-domain simulation, it belongs in the Chapter 20
notebook as supporting synapse material.

### RTM_PLOT_Q

Expose the two-variable synaptic-gate simulation with `q` included in the
state monitor. Use the Python reference parameters and produce the voltage,
`q`, and `s` traces needed for comparison and plotting.

### RTM_PLOT_S_PRESCRIBE_TAU_PEAK

Add deterministic helpers for calculating the synaptic peak time and solving
for `tau_d_q`. Run the two prescribed-peak cases from the Python reference and
return the voltage and synaptic-gate traces. The plotting cell recreates the
three-panel figure.

### RTM_WITH_AUTAPSE_F_I_CURVE

Add a Brian2 RTM neuron with inhibitory autaptic feedback and a callable
current-sweep function. Preserve the Python reference's forward/backward sweep
semantics and returned critical-current values. Keep test parameters
configurable so a focused test can use the smallest sweep that still detects
broken hysteresis or firing-rate behavior, while the plotting cell generates
the full reference figure.

## Validation

Use test-driven development for each example: add a failing Brian-specific
pytest, confirm the expected failure, implement the minimum notebook behavior,
then rerun the focused test.

The Python ports are the primary numeric oracle because they have already been
validated against MATLAB. This avoids repeated MATLAB startup and long MATLAB
sweeps while retaining the existing provenance chain:

`MATLAB reference -> verified Python port -> Brian2 port`

Validation consists of:

- direct array comparison for `B_JAHR_STEVENS`;
- common-grid RMSE comparisons for voltage and synaptic-gate traces;
- comparison of prescribed `tau_d_q` values and resulting peak behavior;
- firing-rate and critical-current comparisons for the autapse sweep, using
  tolerances appropriate for different numerical integration methods;
- regenerated figures visually checked against the existing Python PNGs and
  MATLAB PDFs.

After focused tests pass, run all Chapter 20 Brian tests and then the complete
pytest suite. Update only the four Chapter 20 Brian checklist cells and the
Brian total in `PROGRESS.md`.

## Scope Boundaries

- Do not change the verified Python or MATLAB implementations.
- Do not restructure existing Chapter 20 notebook sections unnecessarily.
- Do not modify or stage existing changes in chapters 19, 33, or 34.
- Do not mark an example complete until both numeric and visual checks pass.
