# Python execution and Brian chapter coverage

The user requested full execution of the 38 Python notebooks and representative,
validated examples for every missing chapter that fits Brian2. They selected
representative coverage, rather than parity with every Python example, and asked
for questions only at the start. Subsequent implementation choices are delegated.

Use three cycles: baseline execution and reference assessment; implementation and
root-cause fixes; complete execution, numerical comparison, and independent review.

Add chapters 10–18, 22, 25–29, 32–33, and 35–40 (23 notebooks). Chapter 03 contains
only static gating functions, so retain its Python implementation and explain the
exclusion in the Brian coverage index. Existing Brian notebooks remain supported.
Each new notebook names its reference, chosen example, units, and omitted portions.
Analytical roots/nullclines remain NumPy/SciPy postprocessing of Brian simulations.

Approaches considered: copy equations into every notebook (easy to inspect but
repeats maintenance); wrap Python reference simulations (does not supply Brian2
implementations); shared Brian2 modules with short teaching notebooks (selected
implementation approach). The user expressly permitted shared modules.

Shared modules live under mnd/brian: neuron factories; low-dimensional dynamics;
phase protocols; common drive/network/plasticity helpers. Use explicit networks
and clocks, local namespaces, numerical millisecond/millivolt conventions explained
at the API boundary, fixed seeds, and finite-state checks. Source kinetic equations,
state initialization, and distinct synaptic gate pools matter more than matching
implementation syntax. Preserve continuous STP/STDP rather than replacing them
with different event-only models. Do not copy the acknowledged chapter27 theta
update bug. Tests may load Python definitions as independent numerical oracles;
production modules and notebooks must not dynamically execute reference notebooks.

Notebooks require this source checkout's editable package; do not pretend the
unreleased helpers are already on PyPI. Keep notebook metadata empty and do not
commit generated outputs. Save executed notebooks, figures, logs, and reports in
ignored results/. Current checkout is retained for the user's IDE and concurrent
validation; changes are divided by file ownership. No new commit or push requested.
