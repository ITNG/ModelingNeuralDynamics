from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]

# No MATLAB-comparison test here: matlab/14/HH_REDUCED_CYCLE_DISTANCE/
# make_figure.m has a variable-name collision (the gating variable `h` is
# later reused as a graphics handle via `h=plot(...)`), which crashes it
# in headless -batch mode partway through panel 2. The checked-in
# figure.pdf must be from an interactive session that tolerated this
# differently. Verified visually instead (matches closely, including the
# panel-to-panel cycle/fixed-point gap shrinking as I decreases).


def _in_window(v, n):
    return (-67.5 < v) & (v < -65.5) & (0.3575 < n) & (n < 0.3675)


def test_hh_reduced_cycle_distance_fixed_points_are_sane():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter14.ipynb")
    panels = ns.simulate_hh_reduced_cycle_distance(t_final=200.0)

    assert [panel[0] for panel in panels] == [5.5, 5.4, 5.3]
    for i_ext, v_c, n_c, v_attr, n_attr, v_rep, n_rep in panels:
        assert -68 < v_c < -65
        assert 0.35 < n_c < 0.38
        # both the settled (attracting) cycle and the reverse-time
        # (repelling-in-forward-time) trajectory actually pass through
        # the zoomed window each panel plots, same as the figure
        assert _in_window(v_attr, n_attr).any()
        assert _in_window(v_rep, n_rep).any()
