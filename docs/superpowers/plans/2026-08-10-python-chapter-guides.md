# Python Chapter Guides Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add practical, self-contained conceptual guides for all 38 existing Python chapter directories, with complete navigation and automated documentation checks.

**Architecture:** Each chapter owns one `README.md` beside its example directories and follows a common reader-oriented template. A central `python/README.md` indexes the guides and learning paths, while `tests/test_python_chapter_docs.py` enforces chapter coverage, required sections, example coverage, and valid relative links.

**Tech Stack:** GitHub-flavored Markdown, Python 3, pytest, `pathlib`, and the standard library.

## Global Constraints

- Use original explanations; do not reproduce the book's prose, exercises, extended derivations, or illustrations.
- Treat the checked-in Python implementation as the source of truth for equations, names, parameter conventions, inputs, and outputs.
- Use the supplied PDF only to verify topic alignment, terminology, and conceptual sequence.
- Target approximately one to two rendered Markdown pages per chapter; chapters with many examples may be longer only to keep the example index complete.
- Use repository-relative links that work locally and on GitHub.
- Document all and only the 38 existing `python/` chapter directories; do not create empty directories for book Chapters 2 or 6.
- Do not change simulations, figures, generated data, example layouts, Brian notebooks, or MATLAB sources.
- Preserve and do not stage unrelated worktree changes.

---

## File Map

- `tests/test_python_chapter_docs.py`: defines the exact documented chapter inventory and validates guide structure, immediate example coverage, and local links.
- `python/README.md`: entry point, complete chapter index, chapter-numbering note, and thematic learning paths.
- `python/<chapter>/README.md`: one practical conceptual companion for the examples in that chapter.
- `README.md`: gains one link to `python/README.md`; installation and Colab content remain unchanged.

The chapter work is split into four reviewer-sized thematic batches. Each batch extends the validator's `CHAPTERS` tuple and adds the matching guides, so every commit is independently complete and testable.

### Task 1: Foundations and Single-Neuron Models (01–09)

**Files:**
- Create: `tests/test_python_chapter_docs.py`
- Create: `python/01_Modeling_a_Single_Neuron/README.md`
- Create: `python/03_The_Classical_HH_ODEs/README.md`
- Create: `python/04_Numerical_Solution_of_HH_ODEs/README.md`
- Create: `python/05_The_Simple_Model_of_Neurons_in_Rodent_Brains/README.md`
- Create: `python/07_Linear_Integrate_and_Fire_(LIF)_Neurons/README.md`
- Create: `python/08_Quadratic_Integrate_and_Fire_(QIF)_and_Theta_Neurons/README.md`
- Create: `python/09_Spike_Frequency_Adaptation/README.md`

**Interfaces:**
- Consumes: immediate example-directory names and model conventions from each chapter's checked-in Python files.
- Produces: the validation contract used by every later documentation batch: `CHAPTERS`, `REQUIRED_HEADINGS`, guide existence/structure checks, example coverage, and relative-link resolution.

- [ ] **Step 1: Add the failing documentation validator**

Create `tests/test_python_chapter_docs.py` with this initial content:

```python
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = REPO_ROOT / "python"
CHAPTERS = (
    "01_Modeling_a_Single_Neuron",
    "03_The_Classical_HH_ODEs",
    "04_Numerical_Solution_of_HH_ODEs",
    "05_The_Simple_Model_of_Neurons_in_Rodent_Brains",
    "07_Linear_Integrate_and_Fire_(LIF)_Neurons",
    "08_Quadratic_Integrate_and_Fire_(QIF)_and_Theta_Neurons",
    "09_Spike_Frequency_Adaptation",
)
REQUIRED_HEADINGS = (
    "## Overview",
    "## Core ideas",
    "## Essential model",
    "## Code examples",
    "## What to look for",
    "## Suggested order",
    "## Prerequisites and related chapters",
    "## Running the examples",
)
LINK_RE = re.compile(r"\[[^\]]+\]\((?:<([^>]+)>|([^)]+))\)")


def guide_path(chapter: str) -> Path:
    return PYTHON_ROOT / chapter / "README.md"


def example_directories(chapter: str) -> list[Path]:
    return sorted(
        path
        for path in (PYTHON_ROOT / chapter).iterdir()
        if path.is_dir()
        and path.name != "__pycache__"
        and not path.name.startswith(".")
    )


@pytest.mark.parametrize("chapter", CHAPTERS)
def test_chapter_guide_has_required_structure(chapter: str) -> None:
    path = guide_path(chapter)
    assert path.is_file(), f"Missing chapter guide: {path.relative_to(REPO_ROOT)}"
    text = path.read_text(encoding="utf-8")
    missing = [heading for heading in REQUIRED_HEADINGS if heading not in text]
    assert not missing, f"{path.relative_to(REPO_ROOT)} missing {missing}"


@pytest.mark.parametrize("chapter", CHAPTERS)
def test_chapter_guide_covers_every_example(chapter: str) -> None:
    path = guide_path(chapter)
    text = path.read_text(encoding="utf-8")
    for example in example_directories(chapter):
        assert f"`{example.name}`" in text, f"Missing description for {example.name}"
        assert f"({example.name}/)" in text, f"Missing link for {example.name}"


@pytest.mark.parametrize("chapter", CHAPTERS)
def test_chapter_guide_relative_links_resolve(chapter: str) -> None:
    path = guide_path(chapter)
    text = path.read_text(encoding="utf-8")
    for match in LINK_RE.finditer(text):
        target = unquote(match.group(1) or match.group(2)).split("#", 1)[0]
        if not target or urlsplit(target).scheme or target.startswith("mailto:"):
            continue
        assert (path.parent / target).resolve().exists(), (
            f"Broken link in {path.relative_to(REPO_ROOT)}: {target}"
        )
```

- [ ] **Step 2: Run the validator and confirm the guides are missing**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: seven failures from `test_chapter_guide_has_required_structure` reporting missing `README.md` files; dependent coverage/link checks may also fail because the files do not yet exist.

- [ ] **Step 3: Inspect the chapter implementations and write the seven guides**

Use the standard headings from `REQUIRED_HEADINGS`; add `## Further reading` only when a useful authoritative reference is available. Define every symbol used in an equation and give one linked bullet per example using the exact form ``- [`EXAMPLE`](EXAMPLE/)``.

Cover these concepts and examples:

- `01`: membrane capacitance, ionic current balance, Hodgkin–Huxley voltage and gating dynamics; `HH_VOLTAGE_TRACE`. Explicitly explain that this repository folder is the first modeling example, not a mirror of the book's vocabulary-only Chapter 1.
- `03`: voltage-dependent activation/inactivation, steady-state gate values, time constants; `HH_GATING_VARIABLES`.
- `04`: numerical integration, transient approach to spiking, limit cycles, and refractoriness; `HH_SOLUTION`, `HH_LIMIT_CYCLE`, `HH_REFRACTORINESS`.
- `05`: the RTM excitatory-neuron model and WB/Erisir inhibitory-neuron models, with their different gating kinetics; `RTM_VOLTAGE_TRACE`, `WB_VOLTAGE_TRACE`, `ERISIR_VOLTAGE_TRACE`, `ERISIR_VOLTAGE_TRACE_2`, `THREE_MODELS_GATING_VARIABLES`.
- `07`: leaky subthreshold integration, threshold/reset, membrane time constant, and comparison with HH; `SUBTHR_FOR_HH`, `TAU_M_FOR_HH`, `LIF_VOLTAGE_TRACE`, `LIF_VOLTAGE_TRACE_2`, `LIF_NEURON_WITH_HH`.
- `08`: QIF finite-time blow-up, reset convention, theta transformation, and circular phase geometry; `QIF_VOLTAGE_TRACE`, `QIF_INFINITE_THRESHOLD`, `THETA_FIRING`, `THREE_CIRCLES`.
- `09`: slow negative feedback, M-current, calcium-dependent AHP current, adaptation maps, and resting/spiking regimes; all nine immediate example directories.

For each chapter, read every example's entry point before describing outputs. Name generated PNG/text files and auxiliary commands only when they are present in that example.

- [ ] **Step 4: Run focused documentation validation**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: `21 passed`.

- [ ] **Step 5: Review prose and commit the batch**

Run: `git diff --check -- tests/test_python_chapter_docs.py python/01_* python/03_* python/04_* python/05_* python/07_* python/08_* python/09_*`

Then stage only the validator and seven guides and commit:

```bash
git add tests/test_python_chapter_docs.py \
  python/01_Modeling_a_Single_Neuron/README.md \
  python/03_The_Classical_HH_ODEs/README.md \
  python/04_Numerical_Solution_of_HH_ODEs/README.md \
  python/05_The_Simple_Model_of_Neurons_in_Rodent_Brains/README.md \
  'python/07_Linear_Integrate_and_Fire_(LIF)_Neurons/README.md' \
  'python/08_Quadratic_Integrate_and_Fire_(QIF)_and_Theta_Neurons/README.md' \
  python/09_Spike_Frequency_Adaptation/README.md
git commit -m "docs: explain foundational neuron examples"
```

### Task 2: Dynamics of Single-Neuron Models (10–19)

**Files:**
- Modify: `tests/test_python_chapter_docs.py`
- Create: `python/10_The_Slow_Fast_Phase_Plane/README.md`
- Create: `python/11_The_Saddle_Node_Bifurcation/README.md`
- Create: `python/12_Two_Dimensional_Bifurcation_Analysis/README.md`
- Create: `python/13_Hopf_Bifurcations/README.md`
- Create: `python/14_Model_Neurons_of_Bifurcation_Type_2/README.md`
- Create: `python/15_Canard_Explosions/README.md`
- Create: `python/16_Model_Neurons_of_Bifurcation_Type_3/README.md`
- Create: `python/17_Frequency_Current_Curves/README.md`
- Create: `python/18_Bistability_Resulting_from_Rebound_Firing/README.md`
- Create: `python/19_Bursting/README.md`

**Interfaces:**
- Consumes: the validator and chapter template established in Task 1.
- Produces: complete documentation for the single-neuron dynamical-systems sequence and expands `CHAPTERS` to 17 entries.

- [ ] **Step 1: Extend the failing chapter inventory**

Append these exact strings to `CHAPTERS`:

```python
    "10_The_Slow_Fast_Phase_Plane",
    "11_The_Saddle_Node_Bifurcation",
    "12_Two_Dimensional_Bifurcation_Analysis",
    "13_Hopf_Bifurcations",
    "14_Model_Neurons_of_Bifurcation_Type_2",
    "15_Canard_Explosions",
    "16_Model_Neurons_of_Bifurcation_Type_3",
    "17_Frequency_Current_Curves",
    "18_Bistability_Resulting_from_Rebound_Firing",
    "19_Bursting",
```

- [ ] **Step 2: Run tests and confirm the ten new guides fail**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: existing chapters pass; failures identify the ten new missing guides.

- [ ] **Step 3: Write guides for Chapters 10–14**

Use the standard template and cover:

- `10`: slow-fast separation, nullclines, reduced HH dynamics, FitzHugh–Nagumo geometry, and changing speed along a cycle; all five example directories.
- `11`: fixed-point creation/destruction at a saddle-node collision and the normal-form diagram; `SADDLE_NODE_BIFURCATION`.
- `12`: state explicitly that the repository's historical folder title corresponds to the book's type-1 onset material; explain two-dimensional RTM reduction, fixed points, invariant cycles, and onset through saddle-node-on-invariant-circle behavior; both example directories.
- `13`: supercritical versus subcritical Hopf bifurcations, stability of equilibria/cycles, phase planes, and bifurcation diagrams; all nine examples.
- `14`: type-2 excitability, finite onset frequency, reduced HH and Erisir phase-plane structure, eigenvalues, and repelling cycles; all seven examples.

- [ ] **Step 4: Write guides for Chapters 15–19**

Use the standard template and cover:

- `15`: canard explosions in slow-fast systems, macro/micro parameter views, FitzHugh–Nagumo and reduced-HH examples; all six examples.
- `16`: type-3 excitability, phasic response, persistent sodium/potassium model, and self-exciting theta-neuron geometry; all five examples.
- `17`: computing firing rate versus applied current, onset frequency, continuous/discontinuous branches, bistability, and model comparisons; all thirteen examples.
- `18`: rebound firing, hysteresis and bistability caused by slow gates, contrasting M- and h-currents across HH, Erisir, and RTM models; all sixteen examples.
- `19`: bursting as alternation between silent and spiking fast-subsystem attractors, slow-current feedback, hysteresis loops, and three-dimensional views; all nine examples.

- [ ] **Step 5: Validate and commit the batch**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: `51 passed`.

Run: `git diff --check -- tests/test_python_chapter_docs.py python/1*/README.md`

Stage `tests/test_python_chapter_docs.py` and only the ten new `README.md` files, then commit:

```bash
git add tests/test_python_chapter_docs.py python/1*/README.md
git commit -m "docs: explain single-neuron dynamics chapters"
```

### Task 3: Communication, Entrainment, and Phase Locking (20–29)

**Files:**
- Modify: `tests/test_python_chapter_docs.py`
- Create: one `README.md` in each chapter directory numbered 20 through 29.

**Interfaces:**
- Consumes: the validator and cross-chapter link conventions from Tasks 1–2.
- Produces: communication/synchronization documentation and expands `CHAPTERS` to 27 entries.

- [ ] **Step 1: Extend `CHAPTERS` with the exact 20–29 directory names**

```python
    "20_Chemical_Synapses",
    "21_Gap_Junctions",
    "22_A_Wilson_Cowan_Model_of_an_Oscillatory_E-I_Network",
    "23_Entrainment_by_Excitatory_Input_Pulses",
    "24_Synchronization_by_Fast_Recurrent_Excitation",
    "25_Phase_Response_Curves_(PRCs)",
    "26_Phase_Locking_of_Two_Oscillators",
    "27_Phase_Locking_with_Delays",
    "28_Weakly_Coupled_Oscillators",
    "29_Stability_of_the_Synchronous_State",
```

- [ ] **Step 2: Run tests and confirm the ten new guides fail**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: the first 17 chapters pass; failures identify Chapters 20–29.

- [ ] **Step 3: Write guides for Chapters 20–24**

Cover every immediate example directory and these concepts:

- `20`: conductance-based synaptic current, gating-variable rise/decay, prescribed peak time, NMDA magnesium block, autapses, and temporal buildup.
- `21`: electrical coupling through voltage differences, synchronization versus subthreshold behavior, and reset-threshold effects in LIF and WB networks.
- `22`: Wilson–Cowan firing-rate variables, E–I nullclines, oscillatory feedback, raster/trace interpretation, and the effect of recurrent excitation.
- `23`: periodic forcing, phase locking, n-to-one entrainment, irregular responses, and return-map analysis for WB and LIF neurons.
- `24`: recurrent excitatory pulses, two-cell/network synchronization, asynchronous and splay initial conditions, and heterogeneity.

- [ ] **Step 4: Write guides for Chapters 25–29**

Cover every immediate example directory and these concepts:

- `25`: finite and infinitesimal phase response curves, phase shifts, type-1/type-2 response, theta-neuron PRCs, RTM/WB pulses, and interaction functions.
- `26`: pulse-coupled phase maps, fixed points and their stability, symmetric/asymmetric abstract PRCs, and RTM comparison.
- `27`: propagation delays in phase maps, two-oscillator and three-oscillator locking, and theta-neuron realization.
- `28`: weak-coupling reduction, interaction/difference functions, identical versus heterogeneous oscillators, and stable phase differences.
- `29`: contraction by an inhibitory pulse, condition numbers/sensitivity, LIF versus RTM behavior, and the theta-neuron river picture.

- [ ] **Step 5: Validate and commit the batch**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: `81 passed`.

Run: `git diff --check -- tests/test_python_chapter_docs.py python/2*/README.md`

Stage `tests/test_python_chapter_docs.py` and only the ten new guides, then commit:

```bash
git add tests/test_python_chapter_docs.py python/2*/README.md
git commit -m "docs: explain neuronal communication and phase locking"
```

### Task 4: Network Rhythms and Plasticity (30–40)

**Files:**
- Modify: `tests/test_python_chapter_docs.py`
- Create: one `README.md` in each chapter directory numbered 30 through 40.

**Interfaces:**
- Consumes: validator/template from earlier tasks and conceptual links to Chapters 9, 20, and 23–29.
- Produces: the final 11 chapter guides and the complete 38-entry `CHAPTERS` inventory.

- [ ] **Step 1: Extend `CHAPTERS` with the exact 30–40 directory names**

```python
    "30_The_PING_Model_of_Gamma_Rhythms",
    "31_ING_Rhythms",
    "32_M_Current_PING_and_Poisson_PING",
    "33_M_Current_PING_and_PINB",
    "34_Nested_Gamma_Theta_Rhythms",
    "35_Periodic_Inhibition",
    "36_F_I_Curves_Pulsed_Excitation",
    "37_Thresholding_in_PING",
    "38_Gamma_Coherence",
    "39_Short-Term_Depression_and_Facilitation",
    "40_Spike_Timing-Dependent_Plasticity(STDP)",
```

- [ ] **Step 2: Run tests and confirm the eleven new guides fail**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: the first 27 chapters pass; failures identify Chapters 30–40.

- [ ] **Step 3: Write guides for Chapters 30–34**

Cover every immediate example directory and these concepts:

- `30`: pyramidal–interneuron gamma (PING), two-cell mechanism, E/I population timing, sparse/random connectivity, drive strength, recurrent inhibition/excitation, LFP and raster outputs.
- `31`: interneuron gamma (ING), inhibitory recovery and synchrony, gap junctions, clustering, abstract pulse coupling, and entrainment of excitatory cells.
- `32`: explain that the repository's historical title covers the book's weak-PING material; stochastic Poisson drive, M-current-mediated deterministic weak PING, phase functions, closeups, and clustering.
- `33`: explain that the repository's historical title covers beta-rhythm material; M-current PING, period skipping, PINB, gap junctions, and cell-assembly timing.
- `34`: nested gamma/theta organization, OLM-cell h- and A-currents, pre-OLM dynamics, theta-modulated excitation/inhibition, and E–I–OLM networks.

- [ ] **Step 4: Write guides for Chapters 35–40**

Cover every immediate example directory and these concepts:

- `35`: periodic inhibitory forcing, oscillatory response windows, and firing-rate changes under inhibition.
- `36`: firing-rate curves under pulsed excitation, square-pulse idealization, and comparison with steady-current f–I curves.
- `37`: threshold/non-reset mechanisms in PING and close-up interpretation of suppression/participation boundaries.
- `38`: coherence of gamma responses to pulses, matched versus mismatched timing, Poisson-PING perturbations, and the meaning of population alignment. Exclude `__pycache__` from examples.
- `39`: dynamic synaptic resources/utilization, depression, facilitation, pulse trains, and consequences in RTM/WB neurons.
- `40`: timing-dependent weight change, Abbott–Song rule, adaptation variable, three-cell PING progression, and PING with evolving synapses.

- [ ] **Step 5: Validate and commit the batch**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: `114 passed`.

Run: `git diff --check -- tests/test_python_chapter_docs.py python/3*/README.md python/40_*/README.md`

Stage `tests/test_python_chapter_docs.py` and only the eleven new guides, then commit:

```bash
git add tests/test_python_chapter_docs.py python/3*/README.md python/40_*/README.md
git commit -m "docs: explain network rhythms and plasticity"
```

### Task 5: Python Index, Root Navigation, and Final Audit

**Files:**
- Modify: `tests/test_python_chapter_docs.py`
- Create: `python/README.md`
- Modify: `README.md`

**Interfaces:**
- Consumes: the complete 38-guide inventory from Tasks 1–4.
- Produces: the public documentation entry point and final repository-wide coverage guarantees.

- [ ] **Step 1: Add failing index and inventory tests**

Append these tests to `tests/test_python_chapter_docs.py`:

```python
def test_inventory_matches_all_python_chapter_directories() -> None:
    actual = {
        path.name
        for path in PYTHON_ROOT.iterdir()
        if path.is_dir() and re.match(r"^\d{2}_", path.name)
    }
    assert set(CHAPTERS) == actual
    assert len(CHAPTERS) == 38


def test_python_index_links_every_chapter() -> None:
    index = PYTHON_ROOT / "README.md"
    assert index.is_file()
    text = index.read_text(encoding="utf-8")
    for chapter in CHAPTERS:
        encoded = chapter.replace("(", "%28").replace(")", "%29")
        assert f"({encoded}/README.md)" in text
    assert "Chapters 2 and 6" in text


def test_root_readme_links_python_guide() -> None:
    text = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    assert "(python/README.md)" in text


def test_all_documentation_relative_links_resolve() -> None:
    paths = [PYTHON_ROOT / "README.md", REPO_ROOT / "README.md"]
    paths.extend(guide_path(chapter) for chapter in CHAPTERS)
    for path in paths:
        text = path.read_text(encoding="utf-8")
        for match in LINK_RE.finditer(text):
            target = unquote(match.group(1) or match.group(2)).split("#", 1)[0]
            if not target or urlsplit(target).scheme or target.startswith("mailto:"):
                continue
            assert (path.parent / target).resolve().exists(), (
                f"Broken link in {path.relative_to(REPO_ROOT)}: {target}"
            )
```

- [ ] **Step 2: Run tests and confirm navigation is missing**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: the index and root-link tests fail because `python/README.md` and its root link do not exist yet; all chapter-guide tests pass.

- [ ] **Step 3: Create the Python chapter index**

Write `python/README.md` with:

- a statement that the directory contains Python programs reproducing and extending computational examples from *An Introduction to Modeling Neuronal Dynamics*;
- a note that repository numbering follows historical code organization, is not always an exact table-of-contents mirror, and Chapters 2 and 6 have no Python example directories;
- all 38 chapter links, using percent-encoded parentheses in link destinations;
- five learning-path groups: foundations and single-neuron models (01–09), single-neuron dynamics (10–19), communication (20–22), entrainment/synchronization/rhythms (23–38), and plasticity (39–40);
- a short “How to use these guides” section directing readers to overview, equations, example descriptions, interpretation, and run instructions.

- [ ] **Step 4: Add root navigation**

Near the root README's “Python codes provided by contributors” text, add one sentence with the exact link destination `python/README.md`, for example:

```markdown
See the [practical Python chapter guides](python/README.md) for the concepts,
equations, example map, and expected results for every implemented chapter.
```

Do not rewrite the existing installation, Colab, introduction, or attribution sections.

- [ ] **Step 5: Run the documentation tests and full test suite**

Run: `python3 -m pytest tests/test_python_chapter_docs.py -q`

Expected: `118 passed`.

Run: `python3 -m pytest -q`

Expected: the entire repository suite passes. If an unrelated pre-existing test fails, record the exact test name/output and verify the documentation-only diff did not touch its inputs before proceeding.

- [ ] **Step 6: Perform the final content audit**

Run:

```bash
git diff --check -- README.md python/README.md \
  ':(glob)python/*/README.md' tests/test_python_chapter_docs.py
find python -mindepth 2 -maxdepth 2 -name README.md | wc -l
grep -R -En 'TBD|TODO|PLACEHOLDER|FIXME' python/*/README.md python/README.md
```

Expected: no whitespace errors; count is `38`; placeholder scan prints nothing.

Manually review Chapters 01, 10, 20, 30, and 40 against their entry points. Confirm each guide accurately represents, respectively, a conductance-based neuron, phase-plane analysis, synaptic dynamics, a network rhythm, and plasticity. Check that no passage reproduces source-book prose and that every equation defines its symbols.

- [ ] **Step 7: Commit navigation and final checks**

Stage only `README.md`, `python/README.md`, and `tests/test_python_chapter_docs.py`, then commit:

```bash
git add README.md python/README.md tests/test_python_chapter_docs.py
git commit -m "docs: add Python chapter guide index"
```

- [ ] **Step 8: Verify the committed result**

Run:

```bash
git status --short
git log -5 --oneline
python3 -m pytest tests/test_python_chapter_docs.py -q
```

Expected: only the user's pre-existing unrelated changes remain in `git status`; the five documentation commits appear in the log; all `118` documentation checks pass.
