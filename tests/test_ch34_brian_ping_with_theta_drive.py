from pathlib import Path

from matlab_ref import load_notebook_as_module

ROOT = Path(__file__).resolve().parents[1]

# python/34_Nested_Gamma_Theta_Rhythms/PING_WITH_THETA_DRIVE draws random
# per-cell drive/phase (no fixed seed in the book's script), so this
# compares aggregate firing rate rather than exact spike times -- a
# 500ms reference run gave ~46 Hz for both E and I populations.
PYTHON_E_RATE_HZ = 923 / 40 / 0.5
PYTHON_I_RATE_HZ = 230 / 10 / 0.5


def test_ping_with_theta_drive_rates_match_python():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter34.ipynb")

    e_rate = len(ns.spm_e_drive.t) / 40 / 1.0
    i_rate = len(ns.spm_i_drive.t) / 10 / 1.0

    assert e_rate > 0 and i_rate > 0, "both populations must fire"
    assert abs(e_rate - PYTHON_E_RATE_HZ) < 15, f"E rate {e_rate:.1f} Hz vs python's {PYTHON_E_RATE_HZ:.1f} Hz"
    assert abs(i_rate - PYTHON_I_RATE_HZ) < 15, f"I rate {i_rate:.1f} Hz vs python's {PYTHON_I_RATE_HZ:.1f} Hz"
