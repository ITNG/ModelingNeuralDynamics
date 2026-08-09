# Design: PyPI publishing for `mnd` + Colab support for `brian/`

## Context

`mnd/` is a small installable package (`mnd/core.py` for shared numpy
helpers, `mnd/brian/` reserved for Brian2-only helpers, currently empty)
already built via setuptools per `pyproject.toml`. Today it's only
installable by cloning the repo and running `pip install -e .` / `uv sync`.

Separately, this repo's README links to a different, standalone package
(`mndynamics`, a separate GitHub repo) built specifically to run the
book's code on Colab/Binder. A past session explicitly decided not to
integrate that package into this repo. This design does the opposite of
that decision for a *different* concern: it makes *this* repo's own `mnd`
package installable from PyPI, and removes the `mndynamics` mention from
the README so the two aren't conflated.

The core motivation is `brian/` chapters specifically: each is a single
notebook per chapter. Opening one directly on Colab today requires the
user to know and `pip install` four or five packages by hand
(`brian2`, `neurodynex3`, etc.) before the first cell runs. Publishing
`mnd` to PyPI with those as declared dependencies turns that into one
line.

`python/` (per-sub-example `.py` scripts, not notebooks) is explicitly
**out of scope** — Colab badges only make sense on notebooks, and
converting ~150 scripts to notebooks is a separate, much larger effort
this design does not cover.

## Scope

1. Publish `mnd` to PyPI as `modelingneuraldynamics` (the name already
   reserved in `pyproject.toml`), via a tag-triggered GitHub Actions
   workflow using PyPI trusted publishing.
2. Add an "Open in Colab" badge and a Colab-only install cell to every
   `brian/chapterNN.ipynb`.
3. Update `README.md`: remove the `mndynamics` mention, add a short
   install section and Colab badges next to the `brian/` chapter list.

Out of scope: converting `python/` to notebooks; populating
`mnd/brian/` with real helpers (still "only when 3+ chapters duplicate
code," per the existing project philosophy); versioning automation
beyond a manual bump + tag.

## Design

### 1. `pyproject.toml` metadata

Add fields PyPI expects that are currently missing:
- `readme = "README.md"`
- `license = "GPL-3.0-or-later"` (matches the existing `LICENSE` file)
- `authors` (repo owner)
- `classifiers` (Python versions, license, intended audience: science/research)
- `urls.Repository` pointing at the GitHub repo

No dependency changes: `mnd`'s existing dependency list (numpy, scipy,
matplotlib, networkx, brian2, neurodynex3, jupyter, ipykernel) already
covers everything a `brian/` notebook needs, so `pip install
modelingneuraldynamics` alone is sufficient for the Colab install cell.

### 2. `.github/workflows/publish.yml`

- Trigger: `push: tags: ["v*"]`.
- Job 1 (`check-version`): parse the tag (strip leading `v`) and the
  `version` field from `pyproject.toml`; fail the workflow if they
  don't match. Cheap guard against publishing with a stale version.
- Job 2 (`build-and-publish`, needs `check-version`): checkout, build
  sdist+wheel (`python -m build`), then `pypa/gh-action-pypi-publish`
  under the `pypi` GitHub environment with `id-token: write` permission
  (no API token secret — matches the trusted publisher already
  registered on PyPI for this repo, workflow file `publish.yml`,
  environment `pypi`).

Release process for a maintainer: bump `version` in `pyproject.toml`,
commit, `git tag vX.Y.Z`, `git push --tags`. No separate GitHub Release
object required (tag push is the trigger).

### 3. Colab support in `brian/*.ipynb`

For each `brian/chapterNN.ipynb`, insert as the first two cells (before
any existing content):

- A markdown cell with the Colab badge:
  ```markdown
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ITNG/ModelingNeuralDynamics/blob/main/brian/chapterNN.ipynb)
  ```
- A code cell:
  ```python
  import sys
  if "google.colab" in sys.modules:
      %pip install -q modelingneuraldynamics
  ```

This runs only on Colab (checked via the `google.colab` module being
present) and is a no-op when run locally in the existing `mnd`
uv/conda env, so it's safe to leave in permanently rather than needing
to be stripped before/after local runs.

Applied uniformly to all existing `brian/chapter*.ipynb` files
(01, 04, 05, 07, 08, 09, 20, plus any added before this ships).

### 4. `README.md` changes

- Remove the paragraph: *"There is also a standalone Python package of
  the book [here](...) to play with codes on Colab and Binder."*
- Add an "Installation" section: `pip install modelingneuraldynamics`,
  one line noting it's the shared helper package used by `brian/`
  notebooks (and some `python/` chapters), with a link to the PyPI
  project page.
- Next to wherever chapters are listed/described, add Colab badges for
  the `brian/` notebooks (same badge markdown as above, linking to each
  chapter file directly).

## Verification

- `python -m build` succeeds locally and `twine check dist/*` passes,
  before ever tagging.
- A first release (e.g. `v0.1.1`) is tagged and pushed; confirm the
  workflow runs green and the package appears on PyPI at
  `pypi.org/project/modelingneuraldynamics`.
- Open one `brian/` notebook's Colab badge link in an actual Colab
  session (manual check, not automatable in CI) and confirm the
  install cell + rest of the notebook runs top-to-bottom.
- `pytest` still passes locally after the `pyproject.toml` metadata
  changes (no behavior change expected, but confirms nothing broke the
  build backend config).

## Decision log

- Publishing happens in *this* repo, not the separate `mndynamics`
  package — user's explicit call, reversing the general spirit (though
  not the letter) of the earlier "don't integrate mndynamics" decision,
  since this is a different package under this repo's own name, not an
  integration of the other one.
- `python/` stays script-only for now; only `brian/` gets Colab badges.
- Trigger is tag push (`v*`), not GitHub Release or manual dispatch.
- Trusted publishing (OIDC), workflow file `publish.yml`, environment
  `pypi` — already configured on the PyPI side by the user.
