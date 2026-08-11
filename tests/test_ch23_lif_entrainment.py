from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "23/LIF_ENTRAINMENT"


@pytest.mark.slow
def test_lif_entrainment_matches_matlab():
    # matlab overwrites v_pre/v_post each loop iteration -- only the final
    # pulse's values survive to the end of the script
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter23.ipynb")
    segments, reset_lines, dotted_lines = ns.simulate_lif_entrainment()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m", ["v_pre", "v_post"])

    tt, v_pre, v_post = dotted_lines[-1]
    assert np.isclose(v_pre, ref["v_pre"], atol=1e-9)
    assert np.isclose(v_post, ref["v_post"], atol=1e-9)
