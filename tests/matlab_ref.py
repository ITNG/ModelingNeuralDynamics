"""Run a MATLAB script headless and pull variables out of its workspace.

Used to check a Python/Brian2 port against the book's original MATLAB code.
"""
import ast
import importlib.util
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

MATLAB = "/home/ziaee/prog/Matlab/R2020a/bin/matlab"


def load_python_port(path):
    """Import a chapter's main.py with cwd set to its own folder.

    Chapter scripts save figures to relative paths (fig.png) and are meant
    to be run as `cd folder && python main.py`. Some don't guard their
    plotting code behind `if __name__ == "__main__"`, so importing them
    with the wrong cwd would silently (re)write a fig.png wherever pytest
    happens to run from.
    """
    path = Path(path).resolve()
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    cwd = os.getcwd()
    try:
        os.chdir(path.parent)
        spec.loader.exec_module(mod)
    finally:
        os.chdir(cwd)
    return mod


def run_matlab_script(script_dir, script_name, varnames, timeout=120):
    """Run script_name (a MATLAB script, not function) inside script_dir and
    return the requested workspace variables as numpy arrays.

    Skips the test (rather than failing) when MATLAB isn't on this machine
    -- it's a local-only license, not available in CI.
    """
    if not (os.path.exists(MATLAB) or shutil.which("matlab")):
        pytest.skip("MATLAB not available on this machine")
    script_dir = Path(script_dir).resolve()
    with tempfile.TemporaryDirectory() as tmp:
        outfile = Path(tmp) / "out.mat"
        var_list = ", ".join(f"'{v}'" for v in varnames)
        cmd = (
            f"cd('{script_dir}'); run('{script_name}'); "
            f"save('{outfile}', {var_list});"
        )
        subprocess.run(
            [MATLAB, "-batch", cmd],
            check=True, timeout=timeout, capture_output=True, text=True,
        )
        data = loadmat(outfile)
    return {v: np.asarray(data[v]).squeeze() for v in varnames}


def trace_rmse(t_ref, v_ref, t_test, v_test):
    """RMSE between two voltage traces, test interpolated onto the ref grid.

    Point-by-point tolerance doesn't work well for spiking traces: a tiny
    timing offset near a spike upstroke (dv/dt very large) can look like a
    huge error at that single sample even though the trace is essentially
    identical. RMSE over the whole trace averages that out.
    """
    v_test_interp = np.interp(t_ref, t_test, v_test)
    return float(np.sqrt(np.mean((v_test_interp - v_ref) ** 2)))


def spike_times(t, v, threshold=0.0):
    """Upward threshold crossings of v(t), as an array of times.

    For repetitively-spiking traces, comparing spike times (rather than
    trace_rmse's whole-trace RMSE) is the right check: MATLAB's fixed-step
    Euler/Heun and Brian2/scipy's adaptive solvers drift apart by a
    fraction of a step on each cycle, so by the last of a dozen spikes the
    traces are offset by ~1ms -- fine on the limit cycle itself, but that
    offset lands squarely on a spike upstroke (huge dv/dt) and blows up a
    point-by-point or whole-trace RMSE even though every spike is present
    and correctly timed to begin with.
    """
    t = np.asarray(t)
    v = np.asarray(v)
    idx = np.where((v[:-1] < threshold) & (v[1:] >= threshold))[0]
    return t[idx]


def load_notebook_as_module(path):
    """Execute a Brian2 chapter notebook's code cells and return the results.

    Mirrors load_python_port for notebook-based chapters (brian/ keeps some
    chapters as a single notebook rather than per-example folders). Runs
    with a non-interactive matplotlib backend and cwd set to the notebook's
    own folder, then returns a namespace with every top-level name the
    notebook defined (functions, StateMonitors, parameters, ...).
    """
    import nbformat
    import matplotlib
    from types import SimpleNamespace

    matplotlib.use("Agg")

    path = Path(path).resolve()
    nb = nbformat.read(path, as_version=4)
    source = "\n".join(
        "".join(cell["source"]) for cell in nb.cells if cell["cell_type"] == "code"
    )
    namespace = {"__name__": path.stem}
    cwd = os.getcwd()
    try:
        os.chdir(path.parent)
        exec(compile(source, str(path), "exec"), namespace)
    finally:
        os.chdir(cwd)
    return SimpleNamespace(**namespace)


def load_notebook_definitions_as_module(path):
    """Load a notebook's imports and definitions without running examples.

    Chapter notebooks commonly combine reusable functions with top-level
    plotting and Brian2 simulations. This focused loader makes definitions
    available to lightweight tests without executing those expensive examples.
    """
    import matplotlib
    import nbformat
    from types import SimpleNamespace

    matplotlib.use("Agg")

    path = Path(path).resolve()
    nb = nbformat.read(path, as_version=4)
    nodes = []
    for cell in nb.cells:
        if cell["cell_type"] != "code":
            continue
        tree = ast.parse("".join(cell["source"]), filename=str(path))
        nodes.extend(
            node
            for node in tree.body
            if isinstance(
                node,
                (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
            )
        )
    source = ast.Module(body=nodes, type_ignores=[])
    namespace = {"__name__": path.stem}
    cwd = os.getcwd()
    try:
        os.chdir(path.parent)
        exec(compile(ast.fix_missing_locations(source), str(path), "exec"), namespace)
    finally:
        os.chdir(cwd)
    return SimpleNamespace(**namespace)
