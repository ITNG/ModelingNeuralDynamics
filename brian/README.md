# Brian2 chapter notebooks

There are **42 notebooks covering 37 chapters**, including 23 new representative
chapter ports. Each new notebook explains the selected experiment and the source
examples it omits. These are representative simulations, not full figure parity
with every Python chapter.

Chapter **03** stays in [Python](../python/chapter03.ipynb): its examples evaluate
static equilibrium gating curves and time constants, with no trajectory to port.
Chapters 02 and 06 have no Python examples in this repository.

## Run from a source checkout

```bash
uv sync --all-groups
uv run jupyter lab
```

Select the checkout’s `.venv` kernel and run every cell from top to bottom.
The new notebooks import the shared modules below from this checkout. Installing
an older PyPI release alone does not provide these new helpers. Existing Colab
links in earlier notebooks do not imply the new source modules are released.

## Coverage

| Chapter | Notebook | Covered example |
|---|---|---|
| 01 | [Open](chapter01.ipynb) | HH voltage trace |
| 04 | [Open](chapter04.ipynb) | HH solution, phase plane and refractoriness |
| 05 | [Open](chapter05.ipynb) | RTM, WB and Erisir neurons |
| 07 | [Open](chapter07.ipynb) | HH subthreshold dynamics and LIF |
| 08 | [Open](chapter08.ipynb) | QIF and theta neurons |
| 09 | [Open](chapter09.ipynb) | M/AHP and LIF adaptation |
| 10 | [Open](chapter10.ipynb) | FHN slow/fast trajectories and reduced HH |
| 11 | [Open](chapter11.ipynb) | Planar saddle-node trajectories and fixed points |
| 12 | [Open](chapter12.ipynb) | Reduced RTM invariant cycle |
| 13 | [Open](chapter13.ipynb) | Supercritical and saturated subcritical Hopf trajectories |
| 14 | [Open](chapter14.ipynb) | Reduced HH rest/limit-cycle geometry |
| 15 | [Open](chapter15.ipynb) | FHN canard transition and step-size comparison |
| 16 | [Open](chapter16.ipynb) | Self-exciting theta neuron with event feedback |
| 17 | [Open](chapter17.ipynb) | Parallel independent-current WB f–I sweep |
| 18 | [Open](chapter18.ipynb) | HH resting/firing bistability |
| 19 | [Open](chapter19.ipynb) | INaP/IK and Erisir bursting |
| 20 | [Open](chapter20.ipynb) | Chemical synaptic gates and buildup |
| 21 | [Open](chapter21.ipynb) | LIF and WB gap junctions |
| 22 | [Open](chapter22.ipynb) | Wilson–Cowan E–I population rates |
| 23 | [Open](chapter23.ipynb) | Periodic excitation and entrainment |
| 24 | [Open](chapter24.ipynb) | RTM recurrent excitation and synchrony |
| 25 | [Open](chapter25.ipynb) | Theta voltage-kick PRC versus closed form |
| 26 | [Open](chapter26.ipynb) | Two pulse-coupled phase oscillators |
| 27 | [Open](chapter27.ipynb) | Delayed pulse-coupled phase oscillators |
| 28 | [Open](chapter28.ipynb) | Weak-coupling phase-difference ODE |
| 29 | [Open](chapter29.ipynb) | LIF first-spike timing after common inhibition |
| 30 | [Open](chapter30/PING_2_cell_with_s_and_q.ipynb) | Two-cell PING and four network variants |
| 31 | [Open](chapter31.ipynb) | ING networks and E-cell entrainment |
| 32 | [Open](chapter32.ipynb) | M-current PING and independent Poisson drive |
| 33 | [Open](chapter33.ipynb) | PINB with distinct inhibitory gate pools |
| 34 | [Open](chapter34.ipynb) | OLM and nested gamma/theta networks |
| 35 | [Open](chapter35.ipynb) | Periodic versus tonic LIF inhibition |
| 36 | [Open](chapter36.ipynb) | RTM constant versus pulsed-drive f–I |
| 37 | [Open](chapter37.ipynb) | RTM recruitment under periodic inhibition |
| 38 | [Open](chapter38.ipynb) | Two-cell gamma coherence with competing drives |
| 39 | [Open](chapter39.ipynb) | Continuous RTM depression and facilitation |
| 40 | [Open](chapter40.ipynb) | Three-cell PING with continuous bounded STDP |

Chapter 30 also contains [PING 1](chapter30/PING_1/PING_1.ipynb),
[PING 2](chapter30/PING_2/PING_2.ipynb), [PING 3](chapter30/PING_3/PING_3.ipynb),
[PING 4](chapter30/PING_4/PING_4.ipynb), and the simpler
[two-cell example](chapter30/PING_2_cell.ipynb).

## Shared simulation modules

| Module | Purpose |
|---|---|
| [`models`](../mnd/brian/models.py) | HH, RTM, WB and Erisir populations; full or instantaneous sodium activation; stable rates at removable singularities |
| [`dynamics`](../mnd/brian/dynamics.py) | Reduced neurons, FHN, Hopf/saddle-node forms, theta feedback and Wilson–Cowan dynamics |
| [`phases`](../mnd/brian/phases.py) | PRCs, pulse events and delays, weak-coupling dynamics, common inhibition |
| [`networks`](../mnd/brian/networks.py) | Continuous q/s gates, E–I connectivity, separate inhibitory pools, independent Poisson streams |
| [`drives`](../mnd/brian/drives.py) | Mean-normalized periodic excitation/inhibition and rate experiments |
| [`plasticity`](../mnd/brian/plasticity.py) | Source continuous depression/facilitation and smooth bounded voltage-driven STDP |

The shared conductance-neuron API uses numeric mV, microamperes, mS and
microfarads; derivatives explicitly divide by Brian’s `ms`. Abstract normal-form
and phase-oscillator time units are mapped to milliseconds for Brian integration
and are labeled as model time, not physiological calibration. Networks and clocks
are explicit; neuron factories do not reset unrelated Brian state.

The chapter 27 example uses independent phase ODEs and queued arrivals; it does
not copy the acknowledged theta-update bug in one legacy Python helper.
Chapters 39–40 preserve the source smooth voltage-triggered plasticity; they do
not substitute event-only resource or spike-pair update rules.

## Verification

```bash
uv run pytest -m 'brian2 and not slow' -n 3 --dist loadfile
```

Tests compare trajectories, spike times, gating, periodic drive and plasticity
against independent SciPy solvers, closed-form formulas and the Python chapters.
Full executed notebooks, figures and reports from this validation are kept in
ignored `results/python-validation/` and `results/brian-expansion/`.
