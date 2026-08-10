from pathlib import Path

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]

# python/24_Synchronization_by_Fast_Recurrent_Excitation/RTM_E_TO_E_NETWORK_2
# prints this locked frequency after a 10s run (too slow to re-run here on
# every test invocation, so it's captured as a reference value).
PYTHON_LOCKED_FREQUENCY_HZ = 27.992475503528404


def test_rtm_e_to_e_network_2_locked_frequency_matches_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter24.ipynb")

    assert abs(ns.frequency - PYTHON_LOCKED_FREQUENCY_HZ) < 1.0, (
        f"locked frequency {ns.frequency:.2f} Hz vs python's {PYTHON_LOCKED_FREQUENCY_HZ:.2f} Hz"
    )
