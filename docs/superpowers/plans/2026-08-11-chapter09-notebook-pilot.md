# Chapter 9 Notebook Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace chapter 9's 8 standalone `main.py` scripts with one Colab-ready notebook (`python/chapter09.ipynb`) with `ipywidgets` sliders, and update its tests to exercise the notebook instead of the scripts.

**Architecture:** Nine HH gating functions duplicated across four of chapter 9's examples move into the `mnd` package. `python/chapter09.ipynb` gets one `simulate_*`/plotting function per example (matching the convention `brian/chapter09.ipynb` already uses) plus an `interact()` widget cell per example that has meaningful tunable parameters. The 6 existing test files switch from importing `main.py` via `load_python_port` to importing function definitions out of the notebook via `load_notebook_definitions_as_module` (already exists in `tests/matlab_ref.py`), so widget cells never execute during test collection.

**Tech Stack:** Python, NumPy, SciPy (`odeint`), Matplotlib, `ipywidgets`, `pytest`, `uv`.

## Global Constraints

- Match existing repo style closely enough to look human-written — same naming/structural conventions as neighboring code, not obviously AI-generated boilerplate.
- Chapter code must stay minimal, simple, self-contained, and beginner-readable (upper-level undergrad audience).
- Numeric behavior of every ported example must be unchanged (verified against the deleted `main.py`'s output, and against MATLAB via the existing `@pytest.mark.slow` tests).
- The notebook's Colab-install cell must use `subprocess.run([sys.executable, "-m", "pip", "install", ...])`, never the `%pip install` IPython magic — `load_notebook_as_module`/`load_notebook_definitions_as_module` `exec(compile(...))` a plain-Python parse of every cell regardless of guards, and `%pip` is a `SyntaxError` outside an IPython kernel (see `scripts/add_colab_cells.py`'s existing comment on this).
- `python/09_Spike_Frequency_Adaptation/README.md` is excluded from `tests/test_python_chapter_docs.py`'s `CHAPTERS` set — do not add it back to that set as part of this plan.

---

## File Map

- `mnd/core.py`: gains 9 HH gating functions (`alpha_h`, `alpha_m`, `alpha_n`, `beta_h`, `beta_m`, `beta_n`, `h_inf`, `m_inf`, `n_inf`), moved verbatim from the four `main.py` files that currently duplicate them.
- `tests/test_mnd_core.py` (new): regression test for the moved gating functions.
- `pyproject.toml`: adds `ipywidgets` as an explicit dependency.
- `python/chapter09.ipynb` (new): the pilot notebook — one function + static figure + (where meaningful) `interact()` widget cell per example.
- `python/09_Spike_Frequency_Adaptation/`: the 8 example folders (`ADAPTATION_MAP`, `CALCIUM_RISE`, `LIF_ADAPT`, `M_CURRENT`, `RTM_AHP`, `RTM_AHP_RESTING`, `RTM_M`, `RTM_M_RESTING`) and their `main.py`/`fig.png` files are deleted.
- `python/09_Spike_Frequency_Adaptation/README.md`: updated to describe the notebook instead of the deleted per-example folders.
- `tests/test_ch09_{rtm_m,rtm_ahp,adaptation_map,calcium_rise,lif_adapt,m_current,v_v_tilde}.py`: rewritten to load `python/chapter09.ipynb` instead of `main.py`.

---

### Task 1: Move HH gating functions into `mnd`

**Files:**
- Modify: `mnd/core.py`
- Test: `tests/test_mnd_core.py`

**Interfaces:**
- Produces: `mnd.core.alpha_h(v)`, `alpha_m(v)`, `alpha_n(v)`, `beta_h(v)`, `beta_m(v)`, `beta_n(v)`, `h_inf(v)`, `m_inf(v)`, `n_inf(v)` — each takes a scalar or NumPy array `v` (membrane voltage in mV) and returns the same shape.

- [ ] **Step 1: Write the failing test**

Create `tests/test_mnd_core.py`:

```python
import pytest

from mnd.core import (
    alpha_h, alpha_m, alpha_n,
    beta_h, beta_m, beta_n,
    h_inf, m_inf, n_inf,
)


def test_gating_functions_at_rest():
    # Reference values at v=-70mV, computed from the original (pre-move)
    # formulas duplicated in python/09_Spike_Frequency_Adaptation/RTM_M
    # and RTM_AHP's main.py.
    v = -70.0
    assert alpha_h(v) == pytest.approx(0.3888296675222378)
    assert alpha_m(v) == pytest.approx(0.09552568506252314)
    assert alpha_n(v) == pytest.approx(0.016180577744981224)
    assert beta_h(v) == pytest.approx(0.000736287619853696)
    assert beta_m(v) == pytest.approx(12.042217041926023)
    assert beta_n(v) == pytest.approx(0.6920153229903757)
    assert h_inf(v) == pytest.approx(0.9981099795551048)
    assert m_inf(v) == pytest.approx(0.00787013592322398)
    assert n_inf(v) == pytest.approx(0.02284760152971809)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_mnd_core.py -v`
Expected: FAIL with `ImportError: cannot import name 'alpha_h' from 'mnd.core'`

- [ ] **Step 3: Add the gating functions to `mnd/core.py`**

Append to `mnd/core.py` (keep the existing `import numpy as np` at the top,
keep the existing `draw_arrow` function unchanged):

```python
def alpha_h(v):
    """Traub-Miles HH-type gating rate."""
    return 0.128 * np.exp(-(v + 50.0) / 18.0)


def alpha_m(v):
    return 0.32 * (v + 54) / (1.0 - np.exp(-(v + 54.0) / 4.0))


def alpha_n(v):
    return 0.032 * (v + 52) / (1.0 - np.exp(-(v + 52.0) / 5.0))


def beta_h(v):
    return 4.0 / (1.0 + np.exp(-(v + 27.0) / 5.0))


def beta_m(v):
    return 0.28 * (v + 27.0) / (np.exp((v + 27.0) / 5.0) - 1.0)


def beta_n(v):
    return 0.5 * np.exp(-(v + 57.0) / 40.0)


def h_inf(v):
    return alpha_h(v) / (alpha_h(v) + beta_h(v))


def m_inf(v):
    return alpha_m(v) / (alpha_m(v) + beta_m(v))


def n_inf(v):
    return alpha_n(v) / (alpha_n(v) + beta_n(v))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/test_mnd_core.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add mnd/core.py tests/test_mnd_core.py
git commit -m "$(cat <<'EOF'
feat: move chapter-9 HH gating functions into mnd.core

alpha_h/alpha_m/alpha_n/beta_h/beta_m/beta_n/h_inf/m_inf/n_inf were
byte-identical across RTM_M, RTM_M_RESTING, RTM_AHP, and
RTM_AHP_RESTING's main.py. First step of the chapter-9 notebook pilot
(see docs/superpowers/specs/2026-08-11-chapter09-notebook-pilot-design.md).
EOF
)"
```

---

### Task 2: Add `ipywidgets` dependency

**Files:**
- Modify: `pyproject.toml`

**Interfaces:**
- Produces: `ipywidgets` importable as a direct (not merely transitive) dependency.

- [ ] **Step 1: Add the dependency**

In `pyproject.toml`, add `"ipywidgets"` to the `dependencies` list (after
`"jupyter"`, before `"ipykernel"`, alphabetically-ish grouping doesn't
matter here — just keep it in that list):

```toml
dependencies = [
    "numpy",
    "scipy",
    "matplotlib",
    "networkx",
    "brian2",
    "jupyter",
    "ipywidgets",
    "ipykernel",
]
```

- [ ] **Step 2: Sync and verify**

Run: `uv sync`
Run: `uv run python -c "import ipywidgets; print(ipywidgets.__version__)"`
Expected: prints a version string (e.g. `8.1.8`), no error.

- [ ] **Step 3: Commit**

```bash
git add pyproject.toml uv.lock
git commit -m "$(cat <<'EOF'
build: add ipywidgets as an explicit dependency

Was only present transitively via jupyter. The chapter-9 notebook
pilot imports it directly for interact() sliders.
EOF
)"
```

---

### Task 3: Build `python/chapter09.ipynb`

**Files:**
- Create: `python/chapter09.ipynb`

**Interfaces:**
- Consumes: `mnd.core.{alpha_h,alpha_m,alpha_n,beta_h,beta_m,beta_n,h_inf,m_inf,n_inf}` (Task 1).
- Produces (top-level names in the notebook, all consumed by Task 4's tests):
  - `m_current(v=None)` → `(v, w_inf_vals, tau_w_vals)`
  - `simulate_rtm_m(i_ext=1.5, g_m=0.25, t_final=300, dt=0.01)` → `(t, v, w)`
  - `calcium_rise(v=None)` → `(v, r)`
  - `simulate_rtm_ahp(i_ext=1.5, g_ahp=0.25, t_final=300, dt=0.01)` → `(t, v, ca)`
  - `simulate_lif_adapt(tau_m=10., I=0.13, tau_a=40., delta=0.05, t_final=300., dt=0.01)` → `(t, v, w)`
  - `adaptation_map(tau_m=10., I=0.12, tau_w=100., delta=0.01, z_max=0.05, dt=0.01, N=100)` → `(z, phi)`
  - `v_v_tilde(tau_m=10., I=0.13, w_k=0.05, tilde_w_k=0.08, tau_w=40., t_final=100., dt=0.01)` → `(t, v, v_tilde)`

This is the largest task in the plan (building the actual notebook content)
and is not further subdivided: the notebook is one deliverable, and a
reviewer evaluates it as a whole (does it run, does it reproduce the old
figures, do the functions have the right signatures) rather than
cell-by-cell.

- [ ] **Step 1: Write the notebook-generation script**

This repo builds notebook JSON directly (see `scripts/add_colab_cells.py`)
rather than via `nbformat`'s object API, so follow that convention. Save
this as a throwaway script (not committed — delete it after Step 3), e.g.
`/tmp/build_chapter09.py`:

```python
import json

def md(*lines):
    return {"cell_type": "markdown", "metadata": {}, "source": [l + "\n" for l in lines[:-1]] + [lines[-1]]}

def code(src):
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": src.splitlines(keepends=True),
    }

cells = []

cells.append(md(
    "# Chapter 9",
    "## Spike Frequency Adaptation",
    "- Code by : [Abolfazl Ziaeemehr](https://github.com/Ziaeemehr)",
))

# ---------------------------------------------------------------- imports
cells.append(code(
    "import numpy as np\n"
    "import matplotlib.pyplot as plt\n"
    "from scipy.integrate import odeint\n"
    "from ipywidgets import interact\n"
    "from mnd.core import alpha_h, alpha_m, alpha_n, beta_h, beta_m, beta_n, h_inf, m_inf, n_inf\n"
))

# ---------------------------------------------------------- M_CURRENT
cells.append(md(
    "## A Model M-Current",
    "",
    "Steady-state activation $w_\\infty(v)$ and time constant $\\tau_w(v)$ "
    "for the slow M-current gate used by the RTM+M-current neuron below.",
))
cells.append(code(
    "def w_inf(v):\n"
    "    return 1.0 / (1.0 + np.exp(-(v + 35.0) / 10.0))\n"
    "\n"
    "\n"
    "def tau_w(v):\n"
    "    return 400.0 / (3.3 * np.exp((v + 35.0) / 20.0) + np.exp(-(v + 35.0) / 20.0))\n"
    "\n"
    "\n"
    "def m_current(v=None):\n"
    "    if v is None:\n"
    "        v = np.arange(-100, 51)\n"
    "    return v, w_inf(v), tau_w(v)\n"
))
cells.append(code(
    "v, wi, tw = m_current()\n"
    "fig, ax = plt.subplots(1, 2, figsize=(10, 4))\n"
    "ax[0].plot(v, wi, color='k', linewidth=2)\n"
    "ax[0].set_xlabel('$v$ [mV]')\n"
    "ax[0].set_ylabel(r'$w_\\infty$')\n"
    "ax[1].plot(v, tw, color='k', linewidth=2)\n"
    "ax[1].set_xlabel('$v$ [mV]')\n"
    "ax[1].set_ylabel(r'$\\tau_w$')\n"
    "plt.tight_layout()\n"
    "plt.show()\n"
))

# ---------------------------------------------------------- RTM_M
cells.append(md(
    "## RTM Neuron with an M-Current",
    "",
    "`m` is not its own state -- the book's `RTM_M/make_figure.m` always "
    "sets `m = m_inf(v)` (quasi-static), never integrates it as an ODE.",
))
cells.append(code(
    "def simulate_rtm_m(i_ext=1.5, g_m=0.25, t_final=300, dt=0.01):\n"
    "    g_k, g_na, g_l = 80, 100, 0.1\n"
    "    v_k, v_na, v_l = -100, 50, -67\n"
    "\n"
    "    def derivative(x0, t):\n"
    "        v, n, h, w = x0\n"
    "        m = m_inf(v)\n"
    "        dv = (i_ext - g_na * h * m ** 3 * (v - v_na) - g_k * n ** 4 * (v - v_k)\n"
    "              - g_l * (v - v_l) - g_m * w * (v - v_k))\n"
    "        dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n\n"
    "        dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h\n"
    "        dw = (w_inf(v) - w) / tau_w(v)\n"
    "        return [dv, dn, dh, dw]\n"
    "\n"
    "    v0 = -70.0\n"
    "    x0 = [v0, n_inf(v0), h_inf(v0), 0.0]\n"
    "    t = np.arange(0, t_final, dt)\n"
    "    sol = odeint(derivative, x0, t)\n"
    "    return t, sol[:, 0], sol[:, 3]\n"
    "\n"
    "\n"
    "def plot_rtm_m(t, v, w):\n"
    "    fig, ax = plt.subplots(2, figsize=(7, 6), sharex=True)\n"
    "    ax[0].plot(t, v, lw=2, c='k')\n"
    "    ax[0].set_ylabel('v [mV]')\n"
    "    ax[1].plot(t, w, lw=2, c='k')\n"
    "    ax[1].set_xlabel('t [ms]')\n"
    "    ax[1].set_ylabel('w')\n"
    "    ax[1].set_ylim(0, max(w) * 1.2)\n"
    "    ax[0].set_xlim(min(t), max(t))\n"
    "    plt.tight_layout()\n"
    "    plt.show()\n"
))
cells.append(code("plot_rtm_m(*simulate_rtm_m())\n"))
cells.append(code(
    "interact(lambda i_ext=1.5, g_m=0.25: plot_rtm_m(*simulate_rtm_m(i_ext=i_ext, g_m=g_m)),\n"
    "         i_ext=(0.0, 3.0, 0.1), g_m=(0.0, 1.0, 0.05));\n"
))
cells.append(md("### RTM Neuron with an M-Current, at Rest ($I_{ext}=0$)"))
cells.append(code("plot_rtm_m(*simulate_rtm_m(i_ext=0, t_final=600))\n"))

# ---------------------------------------------------------- CALCIUM_RISE
cells.append(md(
    "## Calcium-Dependent AHP Currents",
    "",
    "Voltage-dependent calcium target $c_\\infty(v)$ used by the "
    "calcium-dependent AHP current below.",
))
cells.append(code(
    "def cac_inf(v):\n"
    "    return (120 - v) / (1 + np.exp(-(v + 15) / 5)) * 4 / 25\n"
    "\n"
    "\n"
    "def calcium_rise(v=None):\n"
    "    if v is None:\n"
    "        v = np.arange(-1000, 501) / 10.\n"
    "    return v, cac_inf(v)\n"
))
cells.append(code(
    "v, r = calcium_rise()\n"
    "plt.figure(figsize=(7, 4))\n"
    "plt.plot(v, r, color='k', linewidth=2)\n"
    "plt.xlabel('$v$ [mV]')\n"
    "plt.ylabel(r'$c_\\infty$')\n"
    "plt.tight_layout()\n"
    "plt.show()\n"
))

# ---------------------------------------------------------- RTM_AHP
cells.append(md(
    "## RTM Neuron with a Calcium-Dependent AHP Current",
    "",
    "As with the M-current model, `m` is quasi-static "
    "(`m = m_inf(v)`), and the calcium time constant is a plain "
    "constant (80), not voltage-dependent.",
))
cells.append(code(
    "def simulate_rtm_ahp(i_ext=1.5, g_ahp=0.25, t_final=300, dt=0.01):\n"
    "    g_k, g_na, g_l = 80, 100, 0.1\n"
    "    v_k, v_na, v_l = -100, 50, -67\n"
    "\n"
    "    def derivative(x0, t):\n"
    "        v, n, h, ca = x0\n"
    "        m = m_inf(v)\n"
    "        dv = (i_ext - g_na * h * m ** 3 * (v - v_na) - g_k * n ** 4 * (v - v_k)\n"
    "              - g_l * (v - v_l) - g_ahp * ca * (v - v_k))\n"
    "        dn = alpha_n(v) * (1.0 - n) - beta_n(v) * n\n"
    "        dh = alpha_h(v) * (1.0 - h) - beta_h(v) * h\n"
    "        dca = (cac_inf(v) - ca) / 80.0\n"
    "        return [dv, dn, dh, dca]\n"
    "\n"
    "    v0 = -70.0\n"
    "    x0 = [v0, n_inf(v0), h_inf(v0), 0.0]\n"
    "    t = np.arange(0, t_final, dt)\n"
    "    sol = odeint(derivative, x0, t)\n"
    "    return t, sol[:, 0], sol[:, 3]\n"
    "\n"
    "\n"
    "def plot_rtm_ahp(t, v, ca):\n"
    "    fig, ax = plt.subplots(2, figsize=(7, 6), sharex=True)\n"
    "    ax[0].plot(t, v, lw=2, c='k')\n"
    "    ax[0].set_ylabel('v [mV]')\n"
    "    ax[1].plot(t, ca, lw=2, c='k')\n"
    "    ax[1].set_xlabel('t [ms]')\n"
    "    ax[1].set_ylabel('$[Ca^{2+}]$')\n"
    "    ax[1].set_ylim(0, max(ca) * 1.2)\n"
    "    ax[0].set_xlim(min(t), max(t))\n"
    "    plt.tight_layout()\n"
    "    plt.show()\n"
))
cells.append(code("plot_rtm_ahp(*simulate_rtm_ahp())\n"))
cells.append(code(
    "interact(lambda i_ext=1.5, g_ahp=0.25: plot_rtm_ahp(*simulate_rtm_ahp(i_ext=i_ext, g_ahp=g_ahp)),\n"
    "         i_ext=(0.0, 3.0, 0.1), g_ahp=(0.0, 1.0, 0.05));\n"
))
cells.append(md("### RTM Neuron with AHP, at Rest ($I_{ext}=0$)"))
cells.append(code("plot_rtm_ahp(*simulate_rtm_ahp(i_ext=0, t_final=600))\n"))

# ---------------------------------------------------------- LIF_ADAPT
cells.append(md(
    "## LIF Neuron with Spike-Triggered Adaptation",
    "",
    "Dimensionless LIF model (same convention as chapter 7), with an "
    "adaptation current `w` that jumps by `delta` at every spike and "
    "decays with time constant `tau_a`.",
))
cells.append(code(
    "def simulate_lif_adapt(tau_m=10., I=0.13, tau_a=40., delta=0.05, t_final=300., dt=0.01):\n"
    "    def derivative(state):\n"
    "        v, w = state\n"
    "        dv = -v / tau_m + I - w * v\n"
    "        dw = -w / tau_a\n"
    "        return np.array([dv, dw])\n"
    "\n"
    "    def integrate_rk4(x, dt, f):\n"
    "        k1 = dt * f(x)\n"
    "        k2 = dt * f(x + 0.5 * k1)\n"
    "        k3 = dt * f(x + 0.5 * k2)\n"
    "        k4 = dt * f(x + k3)\n"
    "        return x + (k1 + 2. * (k2 + k3) + k4) / 6.\n"
    "\n"
    "    num_steps = int(t_final / dt)\n"
    "    t = np.arange(0, t_final, dt)\n"
    "    v = np.zeros(num_steps)\n"
    "    w = np.zeros(num_steps)\n"
    "\n"
    "    for i in range(1, num_steps):\n"
    "        v_new, w_new = integrate_rk4(np.array([v[i - 1], w[i - 1]]), dt, derivative)\n"
    "        if v_new <= 1:\n"
    "            v[i] = v_new\n"
    "            w[i] = w_new\n"
    "        else:\n"
    "            v[i] = 0.\n"
    "            w[i] = w_new + delta\n"
    "\n"
    "    return t, v, w\n"
    "\n"
    "\n"
    "def plot_lif_adapt(t, v, w):\n"
    "    fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)\n"
    "    ax[0].plot(t, v, color='k', linewidth=2)\n"
    "    ax[0].set_ylabel('$v$')\n"
    "    ax[0].set_ylim(0, 4)\n"
    "    ax[0].set_yticks([0, 1])\n"
    "    ax[1].plot(t, w, color='k', linewidth=2)\n"
    "    ax[1].set_xlabel('$t$')\n"
    "    ax[1].set_ylabel('$w$')\n"
    "    ax[1].set_xlim(0, max(t))\n"
    "    plt.tight_layout()\n"
    "    plt.show()\n"
))
cells.append(code("plot_lif_adapt(*simulate_lif_adapt())\n"))
cells.append(code(
    "interact(lambda I=0.13, tau_a=40., delta=0.05: plot_lif_adapt(*simulate_lif_adapt(I=I, tau_a=tau_a, delta=delta)),\n"
    "         I=(0.05, 0.3, 0.01), tau_a=(10., 100., 5.), delta=(0.0, 0.2, 0.01));\n"
))

# ---------------------------------------------------------- ADAPTATION_MAP
cells.append(md(
    "## Spike-to-Spike Adaptation Map",
    "",
    "$\\phi(z)$: the adaptation value just after the next spike, as a "
    "function of $z$, its value just after the previous spike. This is "
    "an event/threshold-crossing map (Heun/RK2 scheme), not a smooth ODE.",
))
cells.append(code(
    "def adaptation_map(tau_m=10., I=0.12, tau_w=100., delta=0.01, z_max=0.05, dt=0.01, N=100):\n"
    "    dt05 = dt / 2\n"
    "    v = np.zeros(N + 1)\n"
    "    w = np.arange(N + 1) / N * z_max\n"
    "    phi = np.zeros(N + 1)\n"
    "    done = np.zeros(N + 1, dtype=bool)\n"
    "\n"
    "    while not done.all():\n"
    "        v_old = v.copy()\n"
    "        w_old = w.copy()\n"
    "        v_inc = -v / tau_m + I - w * v\n"
    "        w_inc = -w / tau_w\n"
    "        v_tmp = v + dt05 * v_inc\n"
    "        w_tmp = w + dt05 * w_inc\n"
    "        v_inc = -v_tmp / tau_m + I - w_tmp * v_tmp\n"
    "        w_inc = -w_tmp / tau_w\n"
    "        v = v + dt * v_inc\n"
    "        w = w + dt * w_inc\n"
    "\n"
    "        ind = (v > 1) & (~done)\n"
    "        done[ind] = True\n"
    "        phi[ind] = (v[ind] - 1) * w_old[ind] + (1 - v_old[ind]) * w[ind]\n"
    "        phi[ind] = phi[ind] / (v[ind] - v_old[ind]) + delta\n"
    "\n"
    "    z = np.arange(N + 1) / N * z_max\n"
    "    return z, phi\n"
    "\n"
    "\n"
    "def plot_adaptation_map(z, phi, z_max=0.05):\n"
    "    plt.figure(figsize=(6, 6))\n"
    "    plt.plot(z, phi, color='k', linewidth=2)\n"
    "    plt.plot([0, z_max], [0, z_max], color='k', linestyle='dashed')\n"
    "    plt.xlim(0, z_max)\n"
    "    plt.ylim(0, z_max)\n"
    "    plt.gca().set_aspect('equal')\n"
    "    plt.xlabel('$z$')\n"
    "    plt.ylabel(r'$\\phi(z)$')\n"
    "    plt.tight_layout()\n"
    "    plt.show()\n"
))
cells.append(code("plot_adaptation_map(*adaptation_map())\n"))
cells.append(code(
    "interact(lambda I=0.12, delta=0.01: plot_adaptation_map(*adaptation_map(I=I, delta=delta)),\n"
    "         I=(0.05, 0.3, 0.01), delta=(0.0, 0.05, 0.005));\n"
))

# ---------------------------------------------------------- V_V_TILDE
cells.append(md(
    "## $v$ and its Adaptation-Free Counterpart $\\tilde v$",
    "",
    "A purely subthreshold comparison: $v$ has an adaptation-like "
    "multiplicative term $w_k e^{-t/\\tau_w} v$ that decays over time, "
    "shown against a version with a larger initial $\\tilde w_k$. No "
    "threshold or reset here -- both stay subthreshold throughout.",
))
cells.append(code(
    "def v_v_tilde(tau_m=10., I=0.13, w_k=0.05, tilde_w_k=0.08, tau_w=40., t_final=100., dt=0.01):\n"
    "    def derivative(v, t, wk):\n"
    "        return -v / tau_m + I - wk * np.exp(-t / tau_w) * v\n"
    "\n"
    "    t = np.arange(0, t_final + dt, dt)\n"
    "    v = odeint(derivative, 0., t, args=(w_k,))[:, 0]\n"
    "    v_tilde = odeint(derivative, 0., t, args=(tilde_w_k,))[:, 0]\n"
    "    return t, v, v_tilde\n"
    "\n"
    "\n"
    "def plot_v_v_tilde(t, v, v_tilde):\n"
    "    plt.figure(figsize=(7, 4))\n"
    "    plt.plot(t, v, color='k', linewidth=2, label='$v$')\n"
    "    plt.plot(t, v_tilde, color='k', linewidth=2, linestyle='dashed', label=r'$\\tilde{v}$')\n"
    "    plt.xlabel('$t$')\n"
    "    plt.title(r'$v$ (solid) and $\\tilde{v}$ (dashes)')\n"
    "    plt.xlim(0, 60)\n"
    "    plt.tight_layout()\n"
    "    plt.show()\n"
))
cells.append(code("plot_v_v_tilde(*v_v_tilde())\n"))
cells.append(code(
    "interact(lambda w_k=0.05, tilde_w_k=0.08: plot_v_v_tilde(*v_v_tilde(w_k=w_k, tilde_w_k=tilde_w_k)),\n"
    "         w_k=(0.0, 0.2, 0.01), tilde_w_k=(0.0, 0.2, 0.01));\n"
))

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 4,
}

with open("/home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics/python/chapter09.ipynb", "w") as f:
    json.dump(nb, f, indent=1)

print(f"wrote {len(cells)} cells")
```

- [ ] **Step 2: Run the generation script**

Run: `uv run python /tmp/build_chapter09.py`
Expected: `wrote 27 cells`, and `python/chapter09.ipynb` now exists.

- [ ] **Step 3: Add the Colab badge and install cell**

Run: `uv run python scripts/add_colab_cells.py python/chapter09.ipynb`
Expected: `python/chapter09.ipynb: added Colab badge + install cell`

- [ ] **Step 4: Verify the notebook executes end-to-end and matplotlib figures render**

Run: `uv run jupyter nbconvert --to notebook --execute --output /tmp/chapter09_executed.ipynb python/chapter09.ipynb`
Expected: completes with no errors (this actually runs every `simulate_*`
call and every `interact(...)` cell once with its defaults — `interact()`
outside a live kernel just calls the wrapped function once with default
args and doesn't error).

- [ ] **Step 5: Verify numeric outputs match the deleted-to-be `main.py` scripts**

Run this from the repo root (adjust paths if your checkout differs) to
diff every notebook function's output against the still-present `main.py`
scripts before they're deleted in Task 5:

```bash
uv run python3 - <<'PYEOF'
import numpy as np
from scipy.integrate import odeint
from tests.matlab_ref import load_notebook_definitions_as_module, load_python_port
from pathlib import Path

ROOT = Path(".").resolve()
ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")

def old(folder):
    return load_python_port(ROOT / "python" / "09_Spike_Frequency_Adaptation" / folder / "main.py")

o = old("RTM_M")
t, v, w = ns.simulate_rtm_m()
sol = odeint(o.derivative, o.x0, np.arange(0, o.t_final, o.dt))
assert np.allclose(sol[:, 0], v) and np.allclose(sol[:, 3], w), "RTM_M mismatch"

o = old("RTM_M_RESTING")
t, v, w = ns.simulate_rtm_m(i_ext=0, t_final=600)
sol = odeint(o.derivative, o.x0, np.arange(0, o.t_final, o.dt))
assert np.allclose(sol[:, 0], v), "RTM_M_RESTING mismatch"

o = old("RTM_AHP")
t, v, ca = ns.simulate_rtm_ahp()
sol = odeint(o.derivative, o.x0, np.arange(0, o.t_final, o.dt))
assert np.allclose(sol[:, 0], v) and np.allclose(sol[:, 3], ca), "RTM_AHP mismatch"

o = old("RTM_AHP_RESTING")
t, v, ca = ns.simulate_rtm_ahp(i_ext=0, t_final=600)
sol = odeint(o.derivative, o.x0, np.arange(0, o.t_final, o.dt))
assert np.allclose(sol[:, 0], v), "RTM_AHP_RESTING mismatch"

o = old("LIF_ADAPT")
t, v, w = ns.simulate_lif_adapt()
assert np.allclose(o.v, v) and np.allclose(o.w, w), "LIF_ADAPT mismatch"

o = old("ADAPTATION_MAP")
z, phi = ns.adaptation_map()
assert np.allclose(o.phi, phi), "ADAPTATION_MAP mismatch"

o = old("CALCIUM_RISE")
v, r = ns.calcium_rise()
assert np.allclose(o.r, r), "CALCIUM_RISE mismatch"

o = old("M_CURRENT")
v, wi, tw = ns.m_current()
assert np.allclose(o.w_inf, wi) and np.allclose(o.tau_w, tw), "M_CURRENT mismatch"

o = old("V_V_TILDE")
t, v, vt = ns.v_v_tilde()
assert np.allclose(o.v, v) and np.allclose(o.v_tilde, vt), "V_V_TILDE mismatch"

print("all 9 examples match")
PYEOF
```

Expected: `all 9 examples match`.

- [ ] **Step 6: Clean up the throwaway generation script**

```bash
rm /tmp/build_chapter09.py
```

- [ ] **Step 7: Commit**

```bash
git add python/chapter09.ipynb
git commit -m "$(cat <<'EOF'
feat: add python/chapter09.ipynb notebook

One function + static figure + ipywidgets interact() per example,
replacing the 8 main.py scripts (removed in a following commit once
tests are repointed at the notebook). Numeric outputs verified
identical to the scripts they replace.
EOF
)"
```

---

### Task 4: Repoint chapter-9 tests at the notebook

**Files:**
- Modify: `tests/test_ch09_rtm_m.py`
- Modify: `tests/test_ch09_rtm_ahp.py`
- Modify: `tests/test_ch09_adaptation_map.py`
- Modify: `tests/test_ch09_calcium_rise.py`
- Modify: `tests/test_ch09_lif_adapt.py`
- Modify: `tests/test_ch09_m_current.py`
- Modify: `tests/test_ch09_v_v_tilde.py`

**Interfaces:**
- Consumes: `python/chapter09.ipynb`'s functions (Task 3's Produces list) via `tests.matlab_ref.load_notebook_definitions_as_module`.

- [ ] **Step 1: Rewrite `tests/test_ch09_rtm_m.py`**

```python
from pathlib import Path

import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def _check(simulate_kwargs, matlab_dir, v_tol, w_tol):
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    t_p, v_p, w_p = ns.simulate_rtm_m(**simulate_kwargs)

    ref = run_matlab_script(ROOT / "matlab" / matlab_dir, "make_figure.m", ["t", "v", "w"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    rmse_w = trace_rmse(ref["t"], ref["w"], t_p, w_p)
    assert rmse_v < v_tol, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
    assert rmse_w < w_tol, f"w RMSE vs MATLAB too high: {rmse_w:.4f}"


@pytest.mark.slow
def test_rtm_m_matches_matlab():
    # v_tol looser than the resting variant: 8 full-height spikes over
    # 300ms mean even sub-ms timing jitter at the upstrokes contributes
    # non-trivial RMSE, same effect as the ch05 ERISIR trace comparison.
    _check({}, "09/RTM_M", v_tol=15.0, w_tol=0.01)


@pytest.mark.slow
def test_rtm_m_resting_matches_matlab():
    _check(dict(i_ext=0, t_final=600), "09/RTM_M_RESTING", v_tol=5.0, w_tol=0.01)
```

- [ ] **Step 2: Rewrite `tests/test_ch09_rtm_ahp.py`**

```python
from pathlib import Path

import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def _check(simulate_kwargs, matlab_dir, v_tol, ca_tol):
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    t_p, v_p, ca_p = ns.simulate_rtm_ahp(**simulate_kwargs)

    ref = run_matlab_script(ROOT / "matlab" / matlab_dir, "make_figure.m", ["t", "v", "ca"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t_p, v_p)
    rmse_ca = trace_rmse(ref["t"], ref["ca"], t_p, ca_p)
    assert rmse_v < v_tol, f"v RMSE vs MATLAB too high: {rmse_v:.2f}"
    assert rmse_ca < ca_tol, f"ca RMSE vs MATLAB too high: {rmse_ca:.4f}"


@pytest.mark.slow
def test_rtm_ahp_matches_matlab():
    _check({}, "09/RTM_AHP", v_tol=15.0, ca_tol=0.01)


@pytest.mark.slow
def test_rtm_ahp_resting_matches_matlab():
    _check(dict(i_ext=0, t_final=600), "09/RTM_AHP_RESTING", v_tol=5.0, ca_tol=0.001)
```

- [ ] **Step 3: Rewrite `tests/test_ch09_adaptation_map.py`**

```python
from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_adaptation_map_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    z, phi = ns.adaptation_map()

    ref = run_matlab_script(ROOT / "matlab" / "09/ADAPTATION_MAP", "make_figure.m", ["phi"])

    # Same Heun/RK2 scheme and dt as MATLAB (this is an event/threshold-
    # crossing map, not a smooth ODE, so odeint doesn't apply) -- can use a
    # tight-ish tolerance since it's effectively the same discrete algorithm.
    assert np.allclose(phi, ref["phi"], rtol=1e-3, atol=1e-4)
```

- [ ] **Step 4: Rewrite `tests/test_ch09_calcium_rise.py`**

```python
from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_calcium_rise_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    v, r = ns.calcium_rise()

    ref = run_matlab_script(ROOT / "matlab" / "09/CALCIUM_RISE", "make_figure.m", ["r"])

    assert np.allclose(r, ref["r"], rtol=1e-6, atol=1e-6)
```

- [ ] **Step 5: Rewrite `tests/test_ch09_lif_adapt.py`**

```python
from pathlib import Path

import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_lif_adapt_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    t, v, w = ns.simulate_lif_adapt()

    ref = run_matlab_script(ROOT / "matlab" / "09/LIF_ADAPT", "make_figure.m", ["t", "v", "w"])

    rmse_v = trace_rmse(ref["t"], ref["v"], t, v)
    rmse_w = trace_rmse(ref["t"], ref["w"], t, w)
    assert rmse_v < 0.1, f"v RMSE vs MATLAB too high: {rmse_v:.4f}"
    assert rmse_w < 0.01, f"w RMSE vs MATLAB too high: {rmse_w:.4f}"
```

- [ ] **Step 6: Rewrite `tests/test_ch09_m_current.py`**

```python
from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_m_current_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    v, w_inf, tau_w = ns.m_current()

    ref = run_matlab_script(
        ROOT / "matlab" / "09/M_CURRENT", "make_figure.m", ["w_inf", "tau_w"]
    )

    assert np.allclose(w_inf, ref["w_inf"], rtol=1e-6, atol=1e-6)
    assert np.allclose(tau_w, ref["tau_w"], rtol=1e-6, atol=1e-6)
```

- [ ] **Step 7: Rewrite `tests/test_ch09_v_v_tilde.py`**

```python
from pathlib import Path

import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_v_v_tilde_matches_matlab():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter09.ipynb")
    t, v, v_tilde = ns.v_v_tilde()

    # matlab's make_figure.m overwrites its own `v` in place: first with
    # w_k, then again with tilde_w_k. Only the second (tilde_v) survives to
    # be saved -- so that's the only trace matlab and python have both in
    # scope with matching semantics.
    ref = run_matlab_script(ROOT / "matlab" / "09/V_V_TILDE", "make_figure.m", ["v"])

    rmse = trace_rmse(t, ref["v"], t, v_tilde)
    assert rmse < 0.05, f"RMSE vs MATLAB too high: {rmse:.4f}"
```

- [ ] **Step 8: Run the full chapter-9 test suite against MATLAB**

Run: `uv run pytest tests/test_ch09_rtm_m.py tests/test_ch09_rtm_ahp.py tests/test_ch09_adaptation_map.py tests/test_ch09_calcium_rise.py tests/test_ch09_lif_adapt.py tests/test_ch09_m_current.py tests/test_ch09_v_v_tilde.py -v -m ""`

Expected: all 9 tests PASS (7 files, 2 of which have 2 tests each). If
MATLAB isn't installed on this machine, each test SKIPs instead — that's
also acceptable, but if you have a MATLAB license, run this for real
verification before committing, since this is the whole point of the
pilot (numeric correctness against the book's originals).

- [ ] **Step 9: Confirm `interact()` never fires during test collection**

Run: `uv run pytest tests/test_ch09_rtm_m.py -v -m "" -s 2>&1 | grep -i "widget\|comm"`
Expected: no output (no widget-related errors or warnings) — confirms
`load_notebook_definitions_as_module` is correctly skipping the
`interact(...)` cells rather than executing them.

- [ ] **Step 10: Commit**

```bash
git add tests/test_ch09_rtm_m.py tests/test_ch09_rtm_ahp.py tests/test_ch09_adaptation_map.py tests/test_ch09_calcium_rise.py tests/test_ch09_lif_adapt.py tests/test_ch09_m_current.py tests/test_ch09_v_v_tilde.py
git commit -m "$(cat <<'EOF'
test: repoint chapter-9 tests at python/chapter09.ipynb

Switches from load_python_port(main.py) to
load_notebook_definitions_as_module(chapter09.ipynb) -- the latter
strips top-level statements (including the notebook's interact()
widget cells) and keeps only defs/imports, so tests call the
simulate_* functions directly instead of re-running odeint against
module-level derivative/x0.
EOF
)"
```

---

### Task 5: Delete the old per-example folders

**Files:**
- Delete: `python/09_Spike_Frequency_Adaptation/ADAPTATION_MAP/`
- Delete: `python/09_Spike_Frequency_Adaptation/CALCIUM_RISE/`
- Delete: `python/09_Spike_Frequency_Adaptation/LIF_ADAPT/`
- Delete: `python/09_Spike_Frequency_Adaptation/M_CURRENT/`
- Delete: `python/09_Spike_Frequency_Adaptation/RTM_AHP/`
- Delete: `python/09_Spike_Frequency_Adaptation/RTM_AHP_RESTING/`
- Delete: `python/09_Spike_Frequency_Adaptation/RTM_M/`
- Delete: `python/09_Spike_Frequency_Adaptation/RTM_M_RESTING/`
- Modify: `python/09_Spike_Frequency_Adaptation/README.md`

**Interfaces:** none (this task only removes files that Task 4 already
stopped depending on).

- [ ] **Step 1: Confirm nothing else references the folders being deleted**

Run: `grep -rn "09_Spike_Frequency_Adaptation/\(ADAPTATION_MAP\|CALCIUM_RISE\|LIF_ADAPT\|M_CURRENT\|RTM_AHP\|RTM_AHP_RESTING\|RTM_M\|RTM_M_RESTING\)" --include="*.py" --include="*.md" .`

Expected: no matches outside the files you're about to delete/edit (the
grep should be empty, or only match the README.md line you're editing in
Step 3).

- [ ] **Step 2: Delete the folders**

```bash
git rm -r "python/09_Spike_Frequency_Adaptation/ADAPTATION_MAP" \
          "python/09_Spike_Frequency_Adaptation/CALCIUM_RISE" \
          "python/09_Spike_Frequency_Adaptation/LIF_ADAPT" \
          "python/09_Spike_Frequency_Adaptation/M_CURRENT" \
          "python/09_Spike_Frequency_Adaptation/RTM_AHP" \
          "python/09_Spike_Frequency_Adaptation/RTM_AHP_RESTING" \
          "python/09_Spike_Frequency_Adaptation/RTM_M" \
          "python/09_Spike_Frequency_Adaptation/RTM_M_RESTING"
```

- [ ] **Step 3: Update the chapter README**

In `python/09_Spike_Frequency_Adaptation/README.md`, replace the
`## Code examples` section (currently listing the 9 folder links) with:

```markdown
## Code examples

All nine examples now live in one notebook,
[`chapter09.ipynb`](../chapter09.ipynb): `M_CURRENT` and `CALCIUM_RISE`
plot the underlying steady-state laws; `RTM_M`/`RTM_M_RESTING` and
`RTM_AHP`/`RTM_AHP_RESTING` add each slow current to a driven vs.
zero-drive RTM neuron; `LIF_ADAPT` integrates the reset-based adaptation
model; `ADAPTATION_MAP` computes the spike-to-spike map φ(z); and
`V_V_TILDE` compares two subthreshold voltages with different initial
adaptation amplitudes. Each section has an `ipywidgets` slider to explore
its parameters interactively.
```

And replace the `## Running the examples` section with:

```markdown
## Running the examples

Open [`chapter09.ipynb`](../chapter09.ipynb) in Jupyter, or via the Colab
badge at the top of the notebook. Run all cells top to bottom; each
section's static figure reproduces the book's plot, and the `interact(...)`
cell below it lets you adjust that example's parameters with sliders.
```

- [ ] **Step 4: Run the fast test suite to confirm nothing else broke**

Run: `uv run pytest tests/ -v -m "not slow"`
Expected: PASS (this includes `tests/test_python_chapter_docs.py`, which
should still pass since chapter 9's README stays outside its `CHAPTERS`
coverage set).

- [ ] **Step 5: Commit**

```bash
git add python/09_Spike_Frequency_Adaptation/README.md
git commit -m "$(cat <<'EOF'
refactor: remove chapter-9 main.py scripts, superseded by chapter09.ipynb

The 8 per-example folders (main.py + fig.png) are replaced by
python/chapter09.ipynb (added in a prior commit, already the source
tests/test_ch09_*.py exercise). Updates the chapter README to point at
the notebook instead of the deleted folders.
EOF
)"
```

---

### Task 6: Full verification

**Files:** none (verification only).

- [ ] **Step 1: Lint**

Run: `uv run ruff check .`
Expected: no errors.

- [ ] **Step 2: Fast suite**

Run: `uv run pytest tests/ -v -m "not slow"`
Expected: all PASS.

- [ ] **Step 3: Full chapter-9 suite (including MATLAB comparison if available)**

Run: `uv run pytest tests/test_ch09_*.py -v -m ""`
Expected: all PASS (or SKIP for the MATLAB-comparison tests if this
machine has no MATLAB license — but re-run with a license before
considering the pilot done, per this session's earlier finding that
CI itself only runs the full `-m ""` suite on release tags).

- [ ] **Step 4: Confirm the deleted folders are gone and the notebook is the only chapter-9 python source**

Run: `ls python/09_Spike_Frequency_Adaptation/` and `ls python/chapter09.ipynb`
Expected: the folder contains only `README.md`; `python/chapter09.ipynb`
exists.
