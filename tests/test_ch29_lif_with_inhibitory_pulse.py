from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script
import pytest

ROOT = Path(__file__).resolve().parents[1]
MATLAB_DIR = "29/LIF_WITH_INHIBITORY_PULSE"


@pytest.mark.slow
def test_lif_with_inhibitory_pulse_matches_matlab():
    # matlab's script reuses v/period across all N=10 neurons of all 3
    # panels, so by the end of the script they hold only the last
    # neuron's (i=N) trace/period from the third (strong/fast) panel.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter29.ipynb")
    _, _, traces_strong = ns.simulate_lif_pulse_panels()
    ref = run_matlab_script(ROOT / "matlab" / MATLAB_DIR, "make_figure.m",
                             ["v", "period"], timeout=30)

    v_last, k_last, period_last = traces_strong[-1]
    assert np.isclose(period_last, ref["period"], atol=1e-3)
    assert np.allclose(v_last[:k_last + 1], ref["v"][:k_last + 1], atol=1e-3)
