"""Run a MATLAB script headless and pull variables out of its workspace.

Used to check a Python/Brian2 port against the book's original MATLAB code.
"""
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.io import loadmat

MATLAB = "/home/ziaee/prog/Matlab/R2020a/bin/matlab"


def run_matlab_script(script_dir, script_name, varnames, timeout=120):
    """Run script_name (a MATLAB script, not function) inside script_dir and
    return the requested workspace variables as numpy arrays."""
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
