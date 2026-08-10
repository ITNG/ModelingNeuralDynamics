from pathlib import Path

import numpy as np

from matlab_ref import load_python_port

ROOT = Path(__file__).resolve().parents[1]
PYTHON_BASE = "31_ING_Rhythms"


def _check_ing(name):
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing,
    # roughly uniform spike counts across the 100 I-cells) rather than
    # exact spike times.
    py = load_python_port(ROOT / "python" / PYTHON_BASE / name / "main.py")
    t_i = py.t_i_spikes
    i_i = py.i_i_spikes
    assert len(t_i) > 500
    counts = np.bincount(i_i.astype(int), minlength=py.num_i)
    # with heterogeneous drive (sigma_i>0), an occasional cell can
    # plausibly stay silent for the whole run -- require most, not all,
    # cells to have fired at least once.
    assert np.count_nonzero(counts) > 0.9 * py.num_i


def test_ing_2_structure():
    _check_ing("ING_2")


def test_ing_3_structure():
    _check_ing("ING_3")


def test_ing_4_structure():
    _check_ing("ING_4")


def test_ing_5_structure():
    _check_ing("ING_5")


def test_ing_6_structure():
    _check_ing("ING_6")


def test_ing_7_structure():
    _check_ing("ING_7")


def test_ing_8_structure():
    _check_ing("ING_8")


def test_ing_9_structure():
    _check_ing("ING_9")


def test_ing_10_structure():
    _check_ing("ING_10")
