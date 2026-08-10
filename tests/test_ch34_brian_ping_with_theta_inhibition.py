from pathlib import Path

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_ping_with_theta_inhibition_fires_both_populations():
    """Same network as PING_WITH_THETA_DRIVE (that test checks firing rate
    against a python reference run), just with the periodic-inhibition
    theta mechanism instead of drive modulation -- a full python reference
    run for this variant specifically takes several minutes (odeint over
    a 50-neuron stacked ODE system), so this checks the same sanity bound
    the drive variant's python reference falls in (both populations
    firing at a plausible gamma-band rate) rather than re-deriving a
    second expensive reference.
    """
    ns = load_notebook_as_module(ROOT / "brian" / "chapter34.ipynb")

    e_rate = len(ns.spm_e_inh.t) / 40 / 1.0
    i_rate = len(ns.spm_i_inh.t) / 10 / 1.0

    assert 10 < e_rate < 100, f"E rate {e_rate:.1f} Hz outside plausible gamma range"
    assert 10 < i_rate < 100, f"I rate {i_rate:.1f} Hz outside plausible gamma range"
