# PyPI Publishing + Colab Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish this repo's `mnd` package to PyPI as `modelingneuraldynamics` via a tag-triggered GitHub Actions workflow, and add Colab support (badge + one-line install cell) to every `brian/*.ipynb` notebook.

**Architecture:** `pyproject.toml` gains the metadata PyPI requires (readme, license, authors, classifiers, repo URL) with no dependency changes. A new `.github/workflows/publish.yml` builds and publishes `mnd` on `git tag vX.Y.Z` push, using PyPI trusted publishing (OIDC, no stored secret) against an already-configured PyPI trusted publisher (workflow file `publish.yml`, environment `pypi`). Each `brian/*.ipynb` gets two cells prepended (Colab badge markdown + a Colab-only `%pip install` cell). `README.md` drops its mention of the unrelated `mndynamics` package and gains an install section + Colab badges.

**Tech Stack:** Python packaging (`setuptools`, `build`, `twine`), GitHub Actions, `pypa/gh-action-pypi-publish`, Jupyter notebooks (`nbformat`/`jupyter nbconvert`).

## Global Constraints

- PyPI trusted publisher is already registered for this repo with workflow filename `publish.yml` and environment name `pypi` — both must match exactly or publishing fails with an OIDC trust error.
- Publish trigger is a git tag matching `v*` (e.g. `v0.1.1`), not a GitHub Release object and not manual dispatch.
- `python/` is out of scope — no Colab badges or notebook conversion there.
- No new runtime dependencies for `mnd` — `pip install modelingneuraldynamics` must already be sufficient for a `brian/` notebook to run (existing dependency list: numpy, scipy, matplotlib, networkx, brian2, neurodynex3, jupyter, ipykernel).
- Repo remote is `git@github.com:ITNG/ModelingNeuralDynamics.git`; Colab links use `github/ITNG/ModelingNeuralDynamics/blob/main/...`.
- License is GPL-3.0 (per the existing `LICENSE` file, GPLv3 text) — use SPDX identifier `GPL-3.0-or-later`.

---

### Task 1: `pyproject.toml` PyPI metadata + local build verification

**Files:**
- Modify: `pyproject.toml`

**Interfaces:**
- Produces: a `pyproject.toml` `[project]` table that `python -m build` and `twine check` accept without warnings, consumed by Task 2's workflow.

- [ ] **Step 1: Add the missing PyPI metadata fields**

Edit the `[project]` table in `pyproject.toml` from:

```toml
[project]
name = "modelingneuraldynamics"
version = "0.1.0"
description = "Python and Brian2 code companion to 'An Introduction to Modeling Neuronal Dynamics' (Borgers), ported from the book's MATLAB originals"
requires-python = ">=3.10"
dependencies = [
    "numpy",
    "scipy",
    "matplotlib",
    "networkx",
    "brian2",
    "neurodynex3",
    "jupyter",
    "ipykernel",
]
```

to:

```toml
[project]
name = "modelingneuraldynamics"
version = "0.1.0"
description = "Python and Brian2 code companion to 'An Introduction to Modeling Neuronal Dynamics' (Borgers), ported from the book's MATLAB originals"
readme = "README.md"
requires-python = ">=3.10"
license = "GPL-3.0-or-later"
authors = [
    { name = "Abolfazl Ziaeemehr", email = "a.ziaeemehr@gmail.com" },
]
classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering",
]
dependencies = [
    "numpy",
    "scipy",
    "matplotlib",
    "networkx",
    "brian2",
    "neurodynex3",
    "jupyter",
    "ipykernel",
]

[project.urls]
Repository = "https://github.com/ITNG/ModelingNeuralDynamics"
```

Leave every other section of `pyproject.toml` (`[dependency-groups]`, `[build-system]`, `[tool.setuptools.packages.find]`, `[tool.ruff]`, `[tool.ruff.lint]`, `[tool.ruff.lint.per-file-ignores]`) untouched.

- [ ] **Step 2: Build the package locally**

Run: `cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics && /home/ziaee/envs/mnd/bin/python -m pip install --quiet build twine && /home/ziaee/envs/mnd/bin/python -m build`

Expected: `dist/modelingneuraldynamics-0.1.0-py3-none-any.whl` and `dist/modelingneuraldynamics-0.1.0.tar.gz` are created, no errors.

- [ ] **Step 3: Verify the build with twine**

Run: `/home/ziaee/envs/mnd/bin/python -m twine check dist/*`

Expected: `Checking dist/modelingneuraldynamics-0.1.0-py3-none-any.whl: PASSED` and the same for the `.tar.gz`. If it reports warnings about the license format, confirm `license = "GPL-3.0-or-later"` (PEP 639 string form, not the older `{text = "..."}` table form) — that's the form current `setuptools`/`twine` expect.

- [ ] **Step 4: Clean up build artifacts**

Run: `rm -rf dist/ build/ *.egg-info`

These are build-time-only and must not be committed (confirm `git status` shows nothing new before continuing).

- [ ] **Step 5: Run the existing test suite to confirm nothing broke**

Run: `cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics && /home/ziaee/envs/mnd/bin/python -m pytest tests/test_ch01_brian_hh_voltage_trace.py -q`

Expected: `1 passed`. This is a cheap smoke test that the `mnd` package/build backend config change didn't break the existing import machinery used by `tests/matlab_ref.py`.

- [ ] **Step 6: Commit**

```bash
cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics
git add pyproject.toml
git commit -m "Add PyPI metadata to pyproject.toml (readme, license, authors, classifiers)"
```

---

### Task 2: `.github/workflows/publish.yml`

**Files:**
- Create: `.github/workflows/publish.yml`

**Interfaces:**
- Consumes: `pyproject.toml`'s `version` field (Task 1) and the `mnd` package structure already present.
- Produces: nothing consumed by later tasks in this plan — this is the terminal deliverable for the PyPI half of the work.

- [ ] **Step 1: Write the workflow file**

Create `.github/workflows/publish.yml`:

```yaml
name: Publish to PyPI

on:
  push:
    tags:
      - "v*"

jobs:
  check-version:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Verify tag matches pyproject.toml version
        run: |
          TAG_VERSION="${GITHUB_REF_NAME#v}"
          PYPROJECT_VERSION=$(python3 -c "import tomllib; print(tomllib.load(open('pyproject.toml','rb'))['project']['version'])")
          if [ "$TAG_VERSION" != "$PYPROJECT_VERSION" ]; then
            echo "Tag version ($TAG_VERSION) does not match pyproject.toml version ($PYPROJECT_VERSION)"
            exit 1
          fi

  build-and-publish:
    needs: check-version
    runs-on: ubuntu-latest
    environment: pypi
    permissions:
      id-token: write
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - name: Build package
        run: |
          python -m pip install --upgrade build
          python -m build
      - name: Publish to PyPI
        uses: pypa/gh-action-pypi-publish@release/v1
```

Note: `check-version` uses `tomllib`, part of the Python 3.11+ standard library. `ubuntu-latest` GitHub-hosted runners ship Python 3.12 by default as of this writing, so no explicit `setup-python` step is needed in that job — only `build-and-publish` pins Python 3.12 explicitly since it needs a specific interpreter for the build step.

- [ ] **Step 2: Validate the workflow YAML syntax locally**

Run: `cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics && /home/ziaee/envs/mnd/bin/python -c "import yaml; yaml.safe_load(open('.github/workflows/publish.yml'))"`

Expected: no output, no exception (confirms valid YAML). If `yaml` isn't installed in the `mnd` env, run `/home/ziaee/envs/mnd/bin/python -m pip install --quiet pyyaml` first.

- [ ] **Step 3: Commit**

```bash
cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics
git add .github/workflows/publish.yml
git commit -m "Add PyPI publish workflow, triggered on v* tag push"
```

- [ ] **Step 4: Do the first real release (manual, not part of automated task execution)**

This step is a manual action for the repo owner, not something to script:
1. Confirm `pyproject.toml` `version = "0.1.0"` is the version you want to publish first (bump it first if not, repeating Task 1 Step 1-6 with the new version).
2. `git tag v0.1.0 && git push origin v0.1.0`
3. Watch the Actions tab for the `Publish to PyPI` workflow; confirm it goes green.
4. Confirm the package appears at `https://pypi.org/project/modelingneuraldynamics/`.

If this fails with an OIDC/trusted-publisher error, the most likely cause is a mismatch between the workflow filename or environment name here and what's registered on PyPI's trusted publisher form for this project — re-check both against the `publish.yml` name and the `pypi` environment name in the workflow above.

---

### Task 3: Colab badge + install cell in every `brian/*.ipynb`

**Files:**
- Create: `scripts/add_colab_cells.py` (one-off script, not part of the package — used here, can stay in the repo for future new chapters)
- Modify: `brian/chapter01.ipynb`, `brian/chapter04.ipynb`, `brian/chapter05.ipynb`, `brian/chapter07.ipynb`, `brian/chapter08.ipynb`, `brian/chapter09.ipynb`, `brian/chapter20.ipynb`

**Interfaces:**
- Produces: each notebook gains a Colab badge markdown cell + a Colab-install code cell prepended as its new first two cells; the rest of each notebook's cells are unchanged. Consumed by Task 4 (README badges link to these same notebook paths).

- [ ] **Step 1: Write the script that inserts the two cells**

Create `scripts/add_colab_cells.py`:

```python
"""One-off (and reusable-for-new-chapters) script: prepend a Colab badge
and a Colab-only install cell to a brian/*.ipynb notebook.

Usage: python scripts/add_colab_cells.py brian/chapter01.ipynb [more.ipynb ...]
"""
import json
import sys
from pathlib import Path

REPO = "ITNG/ModelingNeuralDynamics"


def badge_cell(notebook_path):
    rel = notebook_path.as_posix()
    url = f"https://colab.research.google.com/github/{REPO}/blob/main/{rel}"
    badge = f"[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({url})"
    return {"cell_type": "markdown", "metadata": {}, "source": [badge]}


def install_cell():
    source = (
        "import sys\n"
        "if \"google.colab\" in sys.modules:\n"
        "    %pip install -q modelingneuraldynamics\n"
    )
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


def add_colab_cells(path):
    path = Path(path)
    nb = json.loads(path.read_text())
    first_source = "".join(nb["cells"][0].get("source", []))
    if "colab.research.google.com" in first_source:
        print(f"{path}: already has a Colab badge, skipping")
        return
    nb["cells"] = [badge_cell(path), install_cell()] + nb["cells"]
    for cell in nb["cells"]:
        cell.pop("id", None)
    path.write_text(json.dumps(nb, indent=1))
    print(f"{path}: added Colab badge + install cell")


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        add_colab_cells(arg)
```

- [ ] **Step 2: Run it against all 7 notebooks**

Run: `cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics && /home/ziaee/envs/mnd/bin/python scripts/add_colab_cells.py brian/chapter01.ipynb brian/chapter04.ipynb brian/chapter05.ipynb brian/chapter07.ipynb brian/chapter08.ipynb brian/chapter09.ipynb brian/chapter20.ipynb`

Expected: 7 lines of `<path>: added Colab badge + install cell`.

- [ ] **Step 3: Verify each notebook still executes cleanly**

Run this for each of the 7 notebooks (substitute the filename):

`/home/ziaee/envs/mnd/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=mnd brian/chapter01.ipynb`

Expected for every file: `[NbConvertApp] Writing ... bytes to brian/chapterNN.ipynb` with no error traceback. The new install cell is a no-op locally (no `google.colab` in `sys.modules` outside actual Colab), so this must produce the exact same simulation outputs as before — only the two new cells at the top are new content.

- [ ] **Step 4: Spot-check the badge URL is well-formed for one notebook**

Run: `/home/ziaee/envs/mnd/bin/python -c "
import json
nb = json.load(open('brian/chapter01.ipynb'))
print(nb['cells'][0]['source'])
print(''.join(nb['cells'][1]['source']))
"`

Expected output:
```
['[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter01.ipynb)']
import sys
if "google.colab" in sys.modules:
    %pip install -q modelingneuraldynamics
```

- [ ] **Step 5: Run the full brian test suite to confirm the cell insertion didn't disturb anything the tests read**

Run: `cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics && /home/ziaee/envs/mnd/bin/python -m pytest tests/ -k "brian" -q`

Expected: all brian-related tests pass (these use `load_notebook_as_module`, which executes every code cell including the two new ones — the install cell's `if` guard means it does nothing outside Colab, so this should behave identically to before).

- [ ] **Step 6: Commit**

```bash
cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics
git add scripts/add_colab_cells.py brian/chapter01.ipynb brian/chapter04.ipynb brian/chapter05.ipynb brian/chapter07.ipynb brian/chapter08.ipynb brian/chapter09.ipynb brian/chapter20.ipynb
git commit -m "Add Colab badge + install cell to all brian/*.ipynb notebooks"
```

---

### Task 4: `README.md` updates

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes: the Colab URLs established in Task 3 (same `blob/main/brian/chapterNN.ipynb` pattern), the PyPI package name from Task 1 (`modelingneuraldynamics`).

- [ ] **Step 1: Remove the `mndynamics` mention**

In `README.md`, delete this paragraph entirely:

```markdown
There is also an standalone Python package of the book [here](https://github.com/Ziaeemehr/mndynamics) to play with codes on Colab and Binder.
```

- [ ] **Step 2: Add an Installation section**

Insert immediately after the paragraph removed in Step 1 (i.e., in the same place, replacing it):

```markdown
### Installation

The shared helper package used by some chapters is on PyPI:

```bash
pip install modelingneuraldynamics
```

### Running chapters on Colab

Every notebook under `brian/` can be opened directly in Google Colab —
click a chapter's badge below. The notebook installs its own
dependencies automatically when running on Colab.

| Chapter | Colab |
|---|---|
| 01 - Modeling a Single Neuron | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter01.ipynb) |
| 04 - Numerical Solution of HH ODEs | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter04.ipynb) |
| 05 - The Simple Model of Neurons in Rodent Brains | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter05.ipynb) |
| 07 - Linear Integrate and Fire (LIF) Neurons | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter07.ipynb) |
| 08 - Quadratic Integrate and Fire (QIF) and Theta Neurons | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter08.ipynb) |
| 09 - Spike Frequency Adaptation | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter09.ipynb) |
| 20 - Chemical Synapses | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapter20.ipynb) |
```

- [ ] **Step 3: Verify no other reference to `mndynamics` remains**

Run: `grep -rn "mndynamics" /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics/README.md`

Expected: no output (no matches).

- [ ] **Step 4: Verify the new table's chapter list matches the notebooks that actually exist**

Run: `ls /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics/brian/*.ipynb`

Expected: exactly `chapter01.ipynb chapter04.ipynb chapter05.ipynb chapter07.ipynb chapter08.ipynb chapter09.ipynb chapter20.ipynb` — matching the 7 rows added in Step 2. If a new `brian/chapterNN.ipynb` was added between writing this plan and executing it, add a corresponding row (same badge pattern) before committing.

- [ ] **Step 5: Commit**

```bash
cd /home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics
git add README.md
git commit -m "README: remove mndynamics mention, add install + Colab badges"
```
