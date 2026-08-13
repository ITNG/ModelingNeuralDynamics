import inspect
from pathlib import Path

import numpy as np
import pytest

from matlab_ref import load_notebook_definitions_as_module, run_matlab_script


ROOT = Path(__file__).resolve().parents[1]


def _ns():
    return load_notebook_definitions_as_module(ROOT / "python" / "chapter38.ipynb")


def test_shared_model_matches_matlab_gate_values():
    ns = _ns()
    v = np.array([-75.0, -60.0, -20.0])

    np.testing.assert_allclose(
        ns.m_e_inf(v), [0.00263048, 0.05624786, 0.94432979], rtol=1e-5
    )
    np.testing.assert_allclose(
        ns.m_i_inf(v), [0.00804324, 0.05293249, 0.81665924], rtol=1e-5
    )
    assert ns.tau_d_q(3.0, 0.3, 0.3) > 0


def test_inhibitory_time_constants_include_matlab_phi_factor():
    ns = _ns()
    v = -60.0
    alpha_h = 0.07 * np.exp(-(v + 58) / 20)
    beta_h = 1 / (np.exp(-0.1 * (v + 28)) + 1)
    alpha_n = -0.01 * (v + 34) / (np.exp(-0.1 * (v + 34)) - 1)
    beta_n = 0.125 * np.exp(-(v + 44) / 80)

    assert ns.tau_h_i(v) == pytest.approx(1 / (alpha_h + beta_h) / 5)
    assert ns.tau_n_i(v) == pytest.approx(1 / (alpha_n + beta_n) / 5)


def test_gamma_coherence_1_produces_aligned_traces_and_spikes():
    ns = _ns()
    result = ns.simulate_gamma_coherence_1(t_final=50.0)

    assert result["t"].shape == result["v_e"].shape
    assert result["t"].shape == result["v_i"].shape
    assert result["t"].shape == result["s_i"].shape
    assert len(result["t_e_spikes"]) > 0
    assert len(result["t_i_spikes"]) > 0


@pytest.mark.slow
def test_gamma_coherence_1_matches_matlab_spike_timing():
    ns = _ns()
    result = ns.simulate_gamma_coherence_1()
    matlab = run_matlab_script(
        ROOT / "matlab" / "38" / "GAMMA_COHERENCE_1",
        "make_figure.m",
        ["t_e_spikes", "t_i_spikes"],
    )

    assert len(result["t_e_spikes"]) == len(matlab["t_e_spikes"])
    assert len(result["t_i_spikes"]) == len(matlab["t_i_spikes"])
    np.testing.assert_allclose(result["t_e_spikes"], matlab["t_e_spikes"], atol=7.0)
    np.testing.assert_allclose(result["t_i_spikes"], matlab["t_i_spikes"], atol=7.0)


def test_gamma_coherence_2_replaces_feedback_with_mean_inhibition():
    ns = _ns()
    result = ns.simulate_gamma_coherence_2(t_final=100.0)

    assert 0.0 < result["mean_s_i"] < 1.0
    assert result["coupled"]["t"].shape == result["mean_inhibition"]["t"].shape
    assert not np.array_equal(
        result["coupled"]["t_e_spikes"],
        result["mean_inhibition"]["t_e_spikes"],
    )


@pytest.mark.slow
def test_gamma_coherence_2_matches_matlab_summary():
    ns = _ns()
    result = ns.simulate_gamma_coherence_2()
    matlab = run_matlab_script(
        ROOT / "matlab" / "38" / "GAMMA_COHERENCE_2",
        "make_figure.m",
        ["t_e_spikes", "t_i_spikes", "mean_s_i"],
    )

    assert abs(result["mean_s_i"] - matlab["mean_s_i"]) < 0.05
    assert abs(len(result["mean_inhibition"]["t_e_spikes"]) - len(matlab["t_e_spikes"])) <= 2
    assert abs(len(result["mean_inhibition"]["t_i_spikes"]) - len(matlab["t_i_spikes"])) <= 2


def test_poisson_ping_is_reproducible_and_active():
    ns = _ns()
    first = ns.simulate_poisson_ping_plus_green(seed=63806, t_final=50.0)
    second = ns.simulate_poisson_ping_plus_green(seed=63806, t_final=50.0)

    bins = np.arange(0, 51, 10)
    first_counts = np.histogram(first["t_e_spikes"], bins)[0]
    second_counts = np.histogram(second["t_e_spikes"], bins)[0]
    assert first["f_hat_e"] == second["f_hat_e"]
    assert first["f_hat_i"] == second["f_hat_i"]
    assert first_counts.mean() == second_counts.mean()
    assert first_counts.std() == second_counts.std()
    assert first["f_hat_e"] > 0
    assert first["f_hat_i"] > 0


def test_poisson_ping_population_statistics_match_matlab():
    ns = _ns()
    result = ns.simulate_poisson_ping_plus_green()
    bins = np.arange(0, 501, 10)
    counts = np.histogram(result["t_e_spikes"], bins)[0]

    # MATLAB R2020a, rng(63806): rates 5.66/31.72 Hz and 10-ms
    # excitatory count mean/std 11.32/13.4245.
    assert abs(result["f_hat_e"] - 5.66) < 1.0
    assert abs(result["f_hat_i"] - 31.72) < 6.0
    assert abs(counts.mean() - 11.32) < 2.0
    assert abs(counts.std() - 13.4245) < 8.0


@pytest.mark.parametrize(
    ("func_name", "period"),
    [
        ("simulate_poisson_ping_plus_pulses", 31.0),
        ("simulate_poisson_ping_mismatched_pulses", 29.0),
    ],
)
def test_poisson_ping_phase_response(func_name, period):
    ns = _ns()
    func = getattr(ns, func_name)
    phases = np.array([0.1, 0.5, 0.9])
    result = func(seed=63806, phases=phases, t_final=30.0)

    assert inspect.signature(func).parameters["period"].default == period
    assert result["spike_counts"].shape == phases.shape
    assert np.all(result["spike_counts"] >= 0)


@pytest.mark.parametrize(
    ("func_name", "matlab_counts"),
    [
        ("simulate_poisson_ping_plus_pulses", np.array([25, 15, 19, 21, 20, 30, 31, 37, 29])),
        ("simulate_poisson_ping_mismatched_pulses", np.array([29, 22, 23, 29, 26, 28, 23, 31, 26])),
    ],
)
def test_poisson_ping_phase_statistics_match_matlab(func_name, matlab_counts):
    ns = _ns()
    result = getattr(ns, func_name)()
    counts = result["spike_counts"]

    assert abs(counts.mean() - matlab_counts.mean()) < 6
    assert abs(counts.std() - matlab_counts.std()) < 4
    assert abs(np.ptp(counts) - np.ptp(matlab_counts)) < 13

    if func_name.endswith("plus_pulses"):
        assert result["fit_slope"] > 0
