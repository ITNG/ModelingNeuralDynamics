from pathlib import Path

import numpy as np
import pytest
from numba import get_num_threads, set_num_threads

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_three_cell_ping_uses_compiled_rhs():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")

    result = ns.simulate_three_cell_ping(ns.reciprocal_g_ee(0.05), t_final=5.0)

    assert result.sol.shape == (250, 15)
    assert np.isfinite(result.sol).all()
    assert ns.derivative_three_cell_ping.nopython_signatures


def test_large_ping_stdp_kernel_parallel_matches_serial_reference():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")
    v_e = np.array([-60.0, -20.0, 10.0])
    g_ee = np.array([[0.01, 0.02, 0.03],
                     [0.04, 0.05, 0.06],
                     [0.07, 0.08, 0.09]])
    a = np.array([1.0, 2.0, 3.0])
    k_plus = np.array([[0.003, 0.004, 0.005],
                       [0.006, 0.007, 0.008],
                       [0.009, 0.010, 0.011]])
    out = np.empty_like(g_ee)

    ns._g_ee_derivative_parallel_s(
        v_e, g_ee, a, 10.0, 10.0, 1.45, k_plus, k_plus * 2 / 3,
        np.full((3, 3), 0.12), np.full((3, 3), 0.005), 3, out,
    )

    expected = np.array([
        [7.4154337100970817e-09, 7.4256667027503234e-05, 4.5469696104530524e-03],
        [-7.4225876738680220e-05, 6.2987658744841307e-05, 1.0415494585136928e-02],
        [-5.4563122585031898e-03, -8.5128062270039932e-03, 5.3897104510543950e-03],
    ])
    np.testing.assert_allclose(out, expected, rtol=1e-13, atol=1e-15)
    assert ns._g_ee_derivative_parallel_s.targetoptions.get("parallel") is True
    overload = next(iter(ns._g_ee_derivative_parallel_s.overloads.values()))
    assert overload.metadata["parfor_diagnostics"].initial_parfors


def test_ping_stdp_compiles_initialization_and_restores_numba_threads():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")
    threads_before = get_num_threads()

    result = ns.simulate_ping_with_stdp(num_e=4, num_i=2, t_final=0.02)
    result_again = ns.simulate_ping_with_stdp(num_e=4, num_i=2, t_final=0.02)

    assert result.g_ee.shape == (4, 4)
    assert np.isfinite(result.g_ee).all()
    np.testing.assert_array_equal(result.g_ee, result_again.g_ee)
    np.testing.assert_array_equal(result.lfp, result_again.lfp)
    np.testing.assert_array_equal(result.t_e_spikes, result_again.t_e_spikes)
    np.testing.assert_array_equal(result.i_e_spikes, result_again.i_e_spikes)
    np.testing.assert_array_equal(result.t_i_spikes, result_again.t_i_spikes)
    np.testing.assert_array_equal(result.i_i_spikes, result_again.i_i_spikes)
    assert ns._rtm_init_population.nopython_signatures
    assert get_num_threads() == threads_before


def test_rtm_init_population_accepts_python_sequences():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")
    i_ext = np.array([0.4, 0.8])
    phi = np.array([0.0, 0.5])

    reference = ns._rtm_init_population.py_func(i_ext, phi)
    result = ns.rtm_init_population(i_ext.tolist(), phi.tolist())

    assert result.shape == (2, 3)
    assert np.isfinite(result).all()
    # njit-compiled and interpreted-Python execution of the same 200,000
    # -step loop can differ by a ULP or two (LLVM vs CPython float codegen),
    # so compare with a tight tolerance rather than bit-for-bit equality.
    np.testing.assert_allclose(result, reference, rtol=0.0, atol=1e-10)


def test_ping_stdp_caps_and_restores_numba_threads_on_error():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")
    function_globals = ns.simulate_ping_with_stdp.__globals__
    original_loop = function_globals["_stdp_loop"]
    threads_before = get_num_threads()
    observed_threads = []

    def raising_loop(*args):
        observed_threads.append(get_num_threads())
        raise RuntimeError("stop after observing the thread mask")

    function_globals["_stdp_loop"] = raising_loop
    try:
        with pytest.raises(RuntimeError, match="thread mask"):
            ns.simulate_ping_with_stdp(num_e=4, num_i=2, t_final=0.02)
        assert get_num_threads() == threads_before
    finally:
        function_globals["_stdp_loop"] = original_loop
        set_num_threads(threads_before)

    assert observed_threads == [min(threads_before, 8)]


def test_three_cell_ping_5_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy (though this particular script draws no random numbers, so
    # results are actually deterministic); verified visually against
    # matlab's figure.pdf (close match: g_12 decays toward 0, g_21
    # grows toward its bound B, and the E1-E2 lag shrinks over time).
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")
    py = ns.simulate_three_cell_ping_5()
    assert len(py.t_e_spikes) > 20
    assert len(py.t_i_spikes) > 10
    assert py.g_12[-1] < py.g_12[0]
    assert py.g_21[-1] > py.g_21[0]
    assert py.g_21[-1] > 0.9 * py.B[1, 0]
    assert len(py.lags) > 5
    assert py.lags[-1] < py.lags[0]


def test_ping_with_stdp_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties rather than exact spikes.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter40.ipynb")
    py = ns.simulate_ping_with_stdp()
    assert len(py.t_e_spikes) > 500
    assert len(py.t_i_spikes) > 100
    assert py.vec_g.min() >= 0.0
    assert py.vec_g.max() <= py.B.max() * 1.01
    # STDP should have produced a spread of synaptic strengths, not
    # left every synapse exactly at its initial value
    assert py.vec_g.std() > 0
