import numpy as np
from matlab_ref import load_notebook_definitions_as_module, REPO_ROOT


def test_periodic_lif_matches_reference_spike_times():
    from mnd.brian.drives import simulate_periodic_lif

    ref = load_notebook_definitions_as_module(REPO_ROOT / "python/chapter35.ipynb")
    _, spikes, rate = ref.run_periodic_inhibition(0.2, True)
    result = simulate_periodic_lif(0.2)
    np.testing.assert_allclose(result.spike_times, spikes, atol=0.03)
    assert result.rate == rate


def test_pinb_has_distinct_pools_and_finite_gates():
    from mnd.brian.networks import simulate_ping

    r = simulate_ping(n_e=4, n_i=2, duration=40.0, tau_ie=20.0)
    assert not np.allclose(r.i_state.s_ie, r.i_state.s_ii)
    for mon in (r.e_state, r.i_state):
        assert np.all(np.isfinite(mon.v))
    assert np.max(r.i_state.s_ie) <= 1


def test_poisson_streams_are_independent_and_reproducible():
    from mnd.brian.networks import simulate_ping

    a = simulate_ping(n_e=4, n_i=2, duration=40.0, poisson_rate=400.0, seed=12)
    b = simulate_ping(n_e=4, n_i=2, duration=40.0, poisson_rate=400.0, seed=12)
    np.testing.assert_array_equal(a.input_spikes.i, b.input_spikes.i)
    np.testing.assert_array_equal(a.input_spikes.t, b.input_spikes.t)
    trains = a.input_spikes.spike_trains()
    assert not np.array_equal(trains[0], trains[1])


def test_pulsed_rtm_rates_match_independent_reference():
    from mnd.brian.drives import simulate_rtm_forcing

    ref = load_notebook_definitions_as_module(REPO_ROOT / "python/chapter36.ipynb")
    currents = np.array([0.3, 0.8, 1.5])
    constant, pulsed, _ = ref.compute_rtm_f_i_curves_pulsed_excitation(
        i_ext_values=currents, t_final=200.0, use_numba=False
    )
    actual_constant = simulate_rtm_forcing(currents, kind="constant", duration=200.0)
    actual_pulsed = simulate_rtm_forcing(currents, kind="pulsed", duration=200.0)
    # At most one crossing can move across an observation boundary.
    np.testing.assert_allclose(actual_constant.rates, constant, atol=5.0)
    np.testing.assert_allclose(actual_pulsed.rates, pulsed, atol=5.0)
    assert np.max(np.abs(actual_constant.v - actual_pulsed.v)) > 20.


def test_two_cell_drive_has_chapter_periods_and_recruits_inhibition():
    from mnd.brian.networks import simulate_two_cell_coherence

    r = simulate_two_cell_coherence(duration=100.0)
    assert len(r.e_spikes.t) > 1 and len(r.i_spikes.t) > 1
    assert np.max(r.i_state.s_ie) > 0.05
    assert np.all(np.isfinite(r.e_state.v))
