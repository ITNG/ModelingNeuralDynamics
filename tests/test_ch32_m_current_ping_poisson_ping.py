from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "32_M_Current_PING_and_Poisson_PING"


def test_plot_phi_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PLOT_PHI" / "main.py")
    assert np.all(py.phi_vec > 0) and np.all(py.phi_vec < 1)
    # phi(x) lies above the diagonal for small x and crosses it once
    assert py.phi_vec[0] > py.x[0]


def test_plot_psi_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PLOT_PSI" / "main.py")
    assert np.all(py.psi_vec > 0) and np.all(py.psi_vec < 1)
    # psi is decreasing (an inhibitory PRC-like map)
    assert py.psi_vec[0] > py.psi_vec[-1]


def test_plot_psi_phi_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PLOT_PSI_PHI" / "main.py")
    assert np.all(py.psi_vec >= 0) and np.all(py.psi_vec <= 1)
    assert np.all(py.phi_vec >= 0) and np.all(py.phi_vec <= 1)


def _check_ping(name, num_e=200, num_i=50, min_e_spikes=100, min_i_spikes=50):
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial E and I
    # firing) rather than exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / name / "main.py")
    assert len(py.t_e_spikes) > min_e_spikes
    assert len(py.t_i_spikes) > min_i_spikes
    assert py.num_e == num_e
    assert py.num_i == num_i


def test_m_current_ping_1_structure():
    py = load_python_port(
        ROOT / "python" / PYTHON_BASE / "M_CURRENT_PING_1" / "main.py"
    )
    # the simulation is an explicit call, not an import side effect
    assert not hasattr(py, "t_e_spikes")
    assert not hasattr(py, "v_plot")
    t_e, i_e, t_i, i_i, v_plot = py.run_smoke(t_final_run=200.0)
    # exact values below hold only for this, the first call to run_smoke()
    # after load_python_port: run_smoke draws from the shared module-level
    # rng, so a second call in this test would silently give different
    # numbers. Do not call run_smoke() again in this test.
    assert len(t_e) == 368
    assert len(t_i) == 282
    assert len(t_e) == len(i_e)
    assert len(t_i) == len(i_i)
    assert np.isclose(v_plot.max(), 48.377, rtol=0.0, atol=0.01)


def test_m_current_ping_1_closeup_structure():
    _check_ping("M_CURRENT_PING_1_CLOSEUP", min_e_spikes=500, min_i_spikes=1000)


def test_m_current_ping_1_from_rest_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "M_CURRENT_PING_1_FROM_REST" / "main.py")
    assert len(py.t_e_spikes) > 500
    assert len(py.t_i_spikes) > 1000
    # network is silent for the first ~200 ms (zero drive, from rest)
    assert py.t_e_spikes.min() > 1.0


def test_m_current_ping_2_closeup_structure():
    _check_ping("M_CURRENT_PING_2_CLOSEUP", min_e_spikes=500, min_i_spikes=1000)


def test_m_current_ping_3_closeup_structure():
    _check_ping("M_CURRENT_PING_3_CLOSEUP", min_e_spikes=500, min_i_spikes=1000)


def test_ping_clusters_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "PING_CLUSTERS" / "main.py")
    assert len(py.t_e_spikes) > 500
    assert len(py.t_i_spikes) > 500
    assert py.num_e == 200
    assert py.num_i == 50


def test_poisson_ping_1_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "POISSON_PING_1" / "main.py")
    assert len(py.t_e_spikes) > 500
    assert len(py.t_i_spikes) > 100
    assert py.num_i == py.num_e // 4


def test_poisson_ping_2_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "POISSON_PING_2" / "main.py")
    assert len(py.t_e_spikes) > 300
    assert len(py.t_i_spikes) > 300


def test_poisson_ping_3_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "POISSON_PING_3" / "main.py")
    assert len(py.t_e_spikes) > 50
    assert len(py.t_i_spikes) > 300
    # strong, all-to-all I-I/E-I coupling with low heterogeneity produces
    # tightly synchronized I-cell volleys
    counts_i = np.bincount(py.i_i_spikes.astype(int), minlength=py.num_i)
    assert np.count_nonzero(counts_i) > 0.9 * py.num_i


def test_poisson_ping_3_voltage_trace_structure():
    py = load_python_port(ROOT / "python" / PYTHON_BASE / "POISSON_PING_3_VOLTAGE_TRACE" / "main.py")
    assert len(py.t_e_spikes) > 50
    assert len(py.t_i_spikes) > 300
    assert len(py.v_one) == py.m_steps + 1
    # the tracked cell spikes at least once (voltage crosses well above rest)
    assert py.v_one.max() > 0
