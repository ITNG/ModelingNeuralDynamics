import pytest

from mnd.core import (
    alpha_h, alpha_m, alpha_n,
    beta_h, beta_m, beta_n,
    h_inf, m_inf, n_inf,
)


def test_gating_functions_at_rest():
    # Reference values at v=-70mV, computed from the original (pre-move)
    # formulas duplicated in python/09_Spike_Frequency_Adaptation/RTM_M
    # and RTM_AHP's main.py.
    v = -70.0
    assert alpha_h(v) == pytest.approx(0.3888296675222378)
    assert alpha_m(v) == pytest.approx(0.09552568506252314)
    assert alpha_n(v) == pytest.approx(0.016180577744981224)
    assert beta_h(v) == pytest.approx(0.000736287619853696)
    assert beta_m(v) == pytest.approx(12.042217041926023)
    assert beta_n(v) == pytest.approx(0.6920153229903757)
    assert h_inf(v) == pytest.approx(0.9981099795551048)
    assert m_inf(v) == pytest.approx(0.00787013592322398)
    assert n_inf(v) == pytest.approx(0.02284760152971809)
