"""Phase protocols checked against analytical and Python chapter references."""
from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module
from mnd.brian.phases import theta_prc, phase_pair, weak_difference, lif_inhibition

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.brian2


def test_theta_prc_matches_closed_form():
    result = theta_prc(phases=np.linspace(.05, .95, 19), dt=.002)
    np.testing.assert_allclose(result.prc, result.predicted_prc, atol=8e-5)
    assert np.isfinite(result.prc).all()
    assert np.max(result.prc) > .1


def test_phase_pair_matches_chapter26_exact_events():
    ref = load_notebook_definitions_as_module(ROOT/'python/chapter26.ipynb')
    a, b, _, _ = ref.simulate_two_pulse_coupled_osc(ref.g_two_pulse_1)
    result = phase_pair(duration=12., dt=.0001)
    for index, expected in enumerate([a[1:], b]):
        expected = expected[expected < 12.-.01]
        actual = result.spike_t[result.spike_i == index][:len(expected)]
        np.testing.assert_allclose(actual, expected, atol=.003)


def test_delayed_phase_pair_matches_chapter27_events():
    ref = load_notebook_definitions_as_module(ROOT/'python/chapter27.ipynb')
    a, b = ref.simulate_two_delayed_pulse_pair(.1, t_final=6.)
    result = phase_pair(kind='delayed', initial_phases=(0., .9), delay=.1,
                        duration=6., dt=.0001, epsilon=1/3, initial_spike=True)
    for index, expected in enumerate([a, b]):
        expected = expected[expected < 6.-.01]
        actual = result.spike_t[result.spike_i == index][:len(expected)]
        np.testing.assert_allclose(actual, expected, atol=.003)


def test_weak_difference_matches_python_trajectory():
    ref = load_notebook_definitions_as_module(ROOT/'python/chapter28.ipynb')
    t, psi = ref.simulate_weakly_coupled_de(ref.wc1_g, .1, .4)
    result = weak_difference(epsilon=.1, initial_phase=.4, duration=60., dt=.01)
    np.testing.assert_allclose(result.psi, np.interp(result.t, t, psi), atol=1e-7)


def test_lif_pulse_first_spikes_match_reference_and_analytic_period():
    ref = load_notebook_definitions_as_module(ROOT/'python/chapter29.ipynb')
    for strength, tau, reference in zip([0., .15, 2.], [1., 9., 1.],
                                      ref.simulate_lif_pulse_panels()):
        result = lif_inhibition(strength=strength, tau_syn=tau)
        expected = np.array([trace[2] for trace in reference])
        np.testing.assert_allclose(result.first_spikes, expected, atol=.002)
        assert np.isfinite(result.v).all()
        if strength == 0:
            np.testing.assert_allclose(result.first_spikes,
                                       (1-result.phases)*result.period, atol=.002)
