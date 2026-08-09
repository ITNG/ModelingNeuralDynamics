# Design: Porting MATLAB chapters to Python/Brian2

## Context

This repo accompanies Christoph Borgers' *An Introduction to Modeling
Neuronal Dynamics*. `matlab/` holds the book's original MATLAB code (38
chapters, 256 individual sub-example scripts). `python/` and `brian/` are
meant to be equivalent ports, for readers who prefer those ecosystems.

A per-chapter-folder audit (comparing sub-example folders, not just chapter
folders) found real completion is much lower than it looked:

- **Python: 59/256 sub-examples (23%)**. Chapters that look "done" at a
  glance are often partial (e.g. ch.17: 3/13, ch.25: 1/13, ch.31: 2/16).
- **Brian2: sparser still.** Existing notebooks (chapters 01, 04, 05, 09,
  20, 30) cover only some sections within their own chapter (e.g.
  `chapter04.ipynb` has 1 of 3 figures, `chapter09.ipynb` has 2 of 9).

The project's core value is staying **minimal and simple enough for
beginners to read** — this governs every decision below.

## Repo/environment setup (already done)

- Reset to `origin`, `master` renamed to `main`, remote repointed to the
  `ITNG` GitHub org (repo was transferred).
- nbstripout wired up via tracked `.gitattributes` (strips notebook
  outputs/exec-counts on commit).
- `.gitignore`: added `.DS_Store`, `__pycache__/`.
- `pyproject.toml`: deps-only (`numpy`, `scipy`, `matplotlib`, `networkx`,
  `brian2`, `neurodynex3`, `jupyter`), `tool.uv.package = false` — no
  installable package, no src-layout. `uv.lock` committed.
- uv env at `/home/ziaee/envs/mnd`, registered as Jupyter kernel `mnd`.
- The separate `mndynamics` package repo is explicitly **not** being
  integrated or referenced.
- This setup is on branch `chore/project-setup`, PR #4 into `main`.

## Code style

- Each MATLAB sub-example is ported as its own self-contained folder
  (`python/NN_ChapterName/EXAMPLE_NAME/main.py[+lib.py]`), mirroring the
  existing convention. Same idea for Brian2, matching whatever structure
  that chapter already uses (single notebook or per-example folder).
- Write it to look human-written: match the naming/structure of
  neighboring chapters, not obviously AI-generated boilerplate.
- Shared code: a real installable package, `mnd/`, flat two-directory
  layout — `mnd/core.py` for plain-numpy helpers shared by `python/`
  chapters (integrators, gating equations), `mnd/brian/` for Brian2-only
  helpers. Chapters `import mnd` after an editable install
  (`uv sync` / `pip install -e .`). Populated incrementally, only with
  code actually duplicated across 3+ chapters — not pre-designed.
  (Superseded the earlier "plain files, no package" call once it became
  clear helper modules needed real organization — see decision log.)
- `pyproject.toml` now builds `mnd` via setuptools
  (`tool.setuptools.packages.find.include = ["mnd*"]`); everything else
  (`python/`, `brian/`, `matlab/`) stays outside the package.
- `brian/input_factory.py` turned out to be a byte-identical vendored
  copy of the already-installed `neurodynex3.tools.input_factory` — not
  moved into `mnd/brian/`, just deleted, and the 5 notebooks that had a
  dead leftover cell importing the local copy (chapters 01, 04, 05, 09,
  20) had that cell removed. Their real import
  (`import neurodynex3.tools.input_factory as input_factory`) already
  worked and was untouched.

## Verification ("definition of done" per sub-example)

Two checks, both required before a sub-example counts as done:

1. **Visual** — regenerate the book's figure(s) for that sub-example.
2. **Numeric** — run the MATLAB original headless (`matlab -batch`,
   MATLAB R2020a is installed locally), dump its key output array(s), and
   add a small `pytest` under `tests/` comparing the Python/Brian2 output
   on a common time grid with a **loose tolerance** (~1e-2 relative).
   Loose because MATLAB uses hand-rolled fixed-step Euler/Heun loops
   while existing Python code already uses `scipy.integrate.odeint`
   (adaptive) — both approximate the same true solution but won't match
   pointwise at tight tolerance. This mismatch predates this effort (e.g.
   `python/04_Numerical_Solution_of_HH_ODEs/HH_SOLUTION/main.py` already
   uses `odeint`), so new chapters keep using `odeint` rather than
   switching to manual Euler to match MATLAB bit-for-bit.

A sub-example is only committed once both checks pass. Never commit
unverified/in-progress work as finished (lesson from the ch.19 /
PING_2-4 work that had to be untracked after being committed prematurely
in PR #4).

## Process

- **Unit of work is the sub-example**, not the chapter (256 total).
- **Order**: strict MATLAB directory order — chapter number ascending,
  then sub-example name within a chapter. Gaps in already-partial
  chapters get filled when we reach that chapter number; no separate
  "finish partials first" pass.
- **Python before Brian2** for a given sub-example — translate the Brian2
  version from the finished Python one when a Python reference exists.
- **Tracking**: a plain `PROGRESS.md` checklist at repo root (matlab item
  → python status → brian status), kept in version control.
- **Branching**: one branch per chapter (or small batch of sub-examples
  within a chapter), PR into `main`.

## Next up

First gap in strict MATLAB order: chapter 05's `ERISIR_VOLTAGE_TRACE_2`
(chapters 01/03/04 are already complete; ch.05 is otherwise done). A
small, low-risk item to validate the workflow (port, visual check,
MATLAB numeric comparison, PROGRESS.md update) before larger chapters
like 09 (9 sub-examples, currently 0 done).
