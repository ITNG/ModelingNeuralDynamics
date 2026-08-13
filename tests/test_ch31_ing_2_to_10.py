from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]
NUM_I = 100


def _check_ing(ns, name):
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties (substantial firing,
    # roughly uniform spike counts across the 100 I-cells) rather than
    # exact spike times.
    t_i, i_i = getattr(ns, f"simulate_{name.lower()}")()
    assert len(t_i) > 500
    counts = np.bincount(i_i.astype(int), minlength=NUM_I)
    # with heterogeneous drive (sigma_i>0), an occasional cell can
    # plausibly stay silent for the whole run -- require most, not all,
    # cells to have fired at least once.
    assert np.count_nonzero(counts) > 0.9 * NUM_I


def _load_ns():
    return load_notebook_definitions_as_module(ROOT / "python" / "chapter31.ipynb")


def test_ing_2_structure():
    _check_ing(_load_ns(), "ING_2")


def test_ing_3_structure():
    _check_ing(_load_ns(), "ING_3")


def test_ing_4_structure():
    _check_ing(_load_ns(), "ING_4")


def test_ing_5_structure():
    _check_ing(_load_ns(), "ING_5")


def test_ing_6_structure():
    _check_ing(_load_ns(), "ING_6")


def test_ing_7_structure():
    _check_ing(_load_ns(), "ING_7")


def test_ing_8_structure():
    _check_ing(_load_ns(), "ING_8")


def test_ing_9_structure():
    _check_ing(_load_ns(), "ING_9")


def test_ing_10_structure():
    _check_ing(_load_ns(), "ING_10")
