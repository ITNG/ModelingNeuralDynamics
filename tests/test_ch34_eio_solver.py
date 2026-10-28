import inspect
from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_eio_excitatory_cell_synapses_follow_reversal_potentials():
    """Each conductance must pull E-cell voltage toward its reversal potential."""
    ns = load_notebook_definitions_as_module(ROOT / "python/chapter34.ipynb")
    # One cell per population; synaptic activation is 0.5 in each population.
    state = np.array([-65., .5, .5, 0, .5] * 2
                     + [-65., .5, .5, .5, .5, .5, 0, .5])
    args = {name: 1. for name in inspect.signature(ns.eio_derivative).parameters
            if name not in ("x0", "t")}
    args.update(
        num_e=1, num_i=1, num_o=1,
        i_ext_e=np.zeros(1), i_ext_i=np.zeros(1), i_ext_o=np.zeros(1),
        v_rev_e=0., v_rev_i=-75., v_rev_o=-75.,
    )
    for name in ("g_ee", "g_ei", "g_eo", "g_ie", "g_ii", "g_io", "g_oe", "g_oi", "g_oo"):
        args[name] = np.zeros((1, 1))
    baseline = ns.eio_derivative(state, 0, **args)
    for coupling, reversal in (("g_ee", 0.), ("g_ie", -75.), ("g_oe", -75.)):
        connected = dict(args, **{coupling: np.ones((1, 1))})
        difference = ns.eio_derivative(state, 0, **connected) - baseline
        np.testing.assert_allclose(difference[0], .5 * (reversal - state[0]),
                                   rtol=0, atol=1e-12)
        np.testing.assert_array_equal(difference[1:], np.zeros(17))
