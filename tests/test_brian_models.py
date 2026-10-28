"""Numerical contracts for the shared Brian neuron factories."""

import brian2 as b2
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from mnd.brian.models import make_neurons
from mnd.core import alpha_h, alpha_n, beta_h, beta_n, h_inf, m_inf, n_inf

pytestmark = pytest.mark.brian2


def test_rtm_trajectory_matches_independent_adaptive_solver():
    initial = [-70.0, h_inf(-70.0), n_inf(-70.0)]

    def rhs(t, y):
        v, h, n = y
        return [
            1.5 + 100 * m_inf(v)**3 * h * (50-v)
            + 80 * n**4 * (-100-v) + .1 * (-67-v),
            alpha_h(v)*(1-h)-beta_h(v)*h,
            alpha_n(v)*(1-n)-beta_n(v)*n,
        ]

    neuron = make_neurons("rtm", current=1.5, dt=.001)
    monitor = b2.StateMonitor(neuron, ["v", "h", "n"], record=True)
    b2.Network(neuron, monitor).run(10*b2.ms)
    t = np.asarray(monitor.t/b2.ms)
    reference = solve_ivp(rhs, (0, 10), initial, t_eval=t, rtol=1e-10, atol=1e-11)
    np.testing.assert_allclose(monitor.v[0], reference.y[0], atol=.002, rtol=0)


@pytest.mark.parametrize("model,voltages", [
    ("hh", [-60., -45.]), ("rtm", [-54., -52., -27.]),
    ("wb", [-35., -34.]), ("erisir", [75.5, 95., -51.25]),
])
def test_rate_singularities_and_short_trajectories_are_finite(model, voltages):
    neuron = make_neurons(model, len(voltages), current=0., v0=voltages, dt=.001)
    monitor = b2.StateMonitor(neuron, ["v", "h", "n"], record=True)
    b2.Network(neuron, monitor).run(.02*b2.ms)
    for field in ["v", "h", "n"]:
        assert np.isfinite(getattr(monitor, field)).all()


def test_factory_has_independent_clock_and_supports_full_sodium_and_extra_current():
    one = make_neurons("rtm", 2, current=[.3, .4], instantaneous_m=False,
                       extra_equations="g_exc : 1", synaptic_current="g_exc*(0-v)",
                       parameters={"C": 2.}, dt=.01)
    two = make_neurons("wb", current=0., dt=.02)
    assert one.clock is not two.clock
    assert one.clock.dt == .01*b2.ms
    assert two.clock.dt == .02*b2.ms
    assert "m" in one.variables
    one.g_exc = .1
    assert np.all(one.I_syn[:] > 0)
    np.testing.assert_allclose(one.I_ext[:], [.3, .4])
    assert one.namespace["C"] == 2.


def test_factory_rejects_invalid_model_and_nonpositive_step():
    with pytest.raises(ValueError, match="model"):
        make_neurons("unknown")
    with pytest.raises(ValueError, match="dt"):
        make_neurons("rtm", dt=0.)
