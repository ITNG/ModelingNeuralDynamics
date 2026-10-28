"""Continuous plasticity is compared against the chapter's independent ODEs."""

import numpy as np
from matlab_ref import load_notebook_definitions_as_module, REPO_ROOT


def test_depression_matches_dynamic_activation_reference():
    from mnd.brian.plasticity import simulate_short_term

    ref = load_notebook_definitions_as_module(REPO_ROOT / "python/chapter39.ipynb")
    # Definition-only loader intentionally excludes the tuple constant assignments.
    ref.derivative_rtm_depressing.__globals__.update(
        g_k=80.0, g_na=100.0, g_l=0.1, v_k=-100.0, v_na=50.0, v_l=-67.0
    )
    expected = ref.simulate_rtm_depressing(t_final=40.0)
    actual = simulate_short_term(duration=40.0)
    for name in ("v", "p", "q", "s"):
        np.testing.assert_allclose(
            getattr(actual, name),
            getattr(expected, name)[: len(actual.t)],
            atol=0.05 if name == "v" else 1e-4,
        )
    assert actual.p.min() < 0.9
    assert np.all((actual.s >= 0) & (actual.s <= 1))


def test_facilitation_changes_release():
    from mnd.brian.plasticity import simulate_short_term

    a = simulate_short_term(duration=100.0, facilitation=True)
    assert a.W.max() > a.W[0] + 0.1
    assert np.all(a.p + a.q <= 1.00001)


def test_stdp_derivative_matches_reference():
    from mnd.brian.plasticity import smooth_weight_derivative

    ref = load_notebook_definitions_as_module(REPO_ROOT / "python/chapter40.ipynb")
    v, a = np.array([20.0, -55.0]), np.array([3.0, 12.0])
    g = np.array([[0.0, 0.05], [0.05, 0.0]])
    out = np.empty_like(g)
    ref._g_ee_derivative_s(
        v, g, a, 10.0, 10.0, 1.45, g, 2 * g / 3, 8 * g, np.maximum(g / 2, 1e-6), 2, out
    )
    np.testing.assert_allclose(smooth_weight_derivative(v, a, g), out, atol=1e-14)


def test_three_cell_weights_evolve_and_stay_bounded():
    from mnd.brian.plasticity import simulate_three_cell_stdp

    a = simulate_three_cell_stdp(duration=100.0)
    assert np.all(np.isfinite(a.weights))
    assert np.all((a.weights >= 0) & (a.weights <= 0.4))
    assert np.max(np.abs(a.weights[:, -1] - 0.05)) > 0.001
    assert len(a.e_spikes.t) and len(a.i_spikes.t)


def test_facilitation_trajectory_matches_python_reference():
    from mnd.brian.plasticity import simulate_short_term

    ref = load_notebook_definitions_as_module(REPO_ROOT / "python/chapter39.ipynb")
    ref.derivative_rtm_depressing_facilitating.__globals__.update(
        g_k=80., g_na=100., g_l=.1, v_k=-100., v_na=50., v_l=-67.
    )
    expected = ref.simulate_rtm_depressing_facilitating(t_final=40.)
    actual = simulate_short_term(duration=40., facilitation=True)
    for name in ("v", "p", "q", "s", "W"):
        np.testing.assert_allclose(
            getattr(actual, name), getattr(expected, name)[:len(actual.t)],
            atol=.05 if name == "v" else 1e-4,
        )


def test_actual_brian_stdp_weight_update_matches_source_at_fixed_states():
    import brian2 as b2
    from mnd.brian.plasticity import simulate_three_cell_stdp

    ref = load_notebook_definitions_as_module(REPO_ROOT / "python/chapter40.ipynb")
    result = simulate_three_cell_stdp(duration=0., dt=1e-5)
    synapses = result.weight_state.source
    # Clamp source/target states; exercise the actual Synapses ODE, including its
    # i=presynaptic/j=postsynaptic matrix convention, without ionic-state drift.
    result.e.active = False
    result.i.active = False
    v, trace = np.array([20., -55.]), np.array([3., 12.])
    result.e.v, result.e.a = v, trace
    weights = np.array([[0., .04], [.08, 0.]])
    synapses.weight = weights[np.asarray(synapses.j), np.asarray(synapses.i)]
    before = np.asarray(synapses.weight[:]).copy()
    initial = np.array([[0., .05], [.05, 0.]])
    expected = np.empty_like(weights)
    ref._g_ee_derivative_s(v, weights, trace, 10., 10., 1.45,
                         initial, 2*initial/3, 8*initial,
                         np.maximum(initial/2, 1e-6), 2, expected)
    result.network.run(1e-5*b2.ms, namespace={})
    actual = (np.asarray(synapses.weight[:])-before)/1e-5
    np.testing.assert_allclose(
        actual, expected[np.asarray(synapses.j), np.asarray(synapses.i)], atol=1e-7,
    )
