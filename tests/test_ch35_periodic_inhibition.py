from pathlib import Path

import numpy as np

from matlab_ref import load_notebook_definitions_as_module

ROOT = Path(__file__).resolve().parents[1]


def test_oscillations_structure():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter35.ipynb")
    phi0 = ns.phi_of(1e-5)
    phi10 = ns.phi_of(10)
    # larger alpha sharpens the pulses: same mean (normalized to 1 in
    # the sense that mean over one period is fixed), but higher peak
    assert phi10.max() > phi0.max()


def test_periodic_inhibition_structure():
    # deterministic script (no RNG) -- exact match confirmed visually
    # against matlab's figure.pdf, so check the exact spike frequencies.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter35.ipynb")
    _, _, freq_015 = ns.run_periodic_inhibition(0.15, True, alpha=5.0)
    _, _, freq_02 = ns.run_periodic_inhibition(0.2, True, alpha=5.0)
    assert freq_015 == 40.0
    assert freq_02 == 80.0
    # tonic (mean) inhibition never lets I=0.15 or 0.2 reach threshold
    _, _, freq_bar_015 = ns.run_periodic_inhibition(0.15, False, alpha=5.0)
    assert freq_bar_015 == 0.0


def test_periodic_inhibition_3_structure():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter35.ipynb")
    _, _, freq_015 = ns.run_periodic_inhibition(0.15, True, alpha=1.0)
    _, _, freq_02 = ns.run_periodic_inhibition(0.2, True, alpha=1.0)
    assert freq_015 == 40.0
    assert freq_02 == 80.0


def test_periodic_inhibition_2_structure():
    # matlab's rng('default'); rng(63806) can't be bit-reproduced by
    # numpy, so we check structural properties instead of exact spikes.
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter35.ipynb")
    _, _, freq = ns.run_periodic_inhibition_noisy(True)
    _, _, freq_bar = ns.run_periodic_inhibition_noisy(False)
    assert freq > freq_bar


def test_periodic_inhibition_f_i_curve_structure():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter35.ipynb")
    # the dense scan is an explicit call, not an import side effect
    assert not hasattr(ns, "f_vec")
    I_vec, f_vec = ns.compute_periodic_inhibition_f_i_curve()
    # visually confirmed step values (40, 80, 120, 160 Hz plateaus)
    assert set(np.round(np.unique(f_vec))) >= {0.0, 40.0, 80.0, 120.0, 160.0}
    assert f_vec[0] == 0.0
    assert f_vec[-1] == 160.0
    # sweep starts at 0.11 (I_vec = arange(1,101)/100*0.2+0.1)
    assert np.isclose(I_vec[0], 0.102)


def test_periodic_inhibition_f_i_curve_2_structure():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter35.ipynb")
    _, f_vec = ns.compute_periodic_inhibition_f_i_curve_2()
    assert set(np.round(np.unique(f_vec))) >= {0.0, 40.0, 80.0, 120.0, 160.0}


def _load_rtm_notebook():
    ns = load_notebook_definitions_as_module(ROOT / "python" / "chapter35.ipynb")
    # the dense scans are explicit calls, not import side effects
    assert not hasattr(ns, "f_vec_tonic")
    assert not hasattr(ns, "f_vec_periodic")
    return ns


def test_rtm_f_i_curve_with_inhibition_structure():
    ns = _load_rtm_notebook()
    i_ext_vec, f_vec_tonic, f_vec_periodic = ns.compute_rtm_f_i_curves()
    assert len(f_vec_tonic) == len(i_ext_vec)
    assert len(f_vec_periodic) == len(i_ext_vec)
    # both curves are silent at low drive and fire at high drive
    assert f_vec_tonic[0] == 0.0
    assert f_vec_tonic[-1] > 0.0
    assert f_vec_periodic[-1] > 0.0


def test_rtm_f_i_curve_with_inhibition_2_structure():
    ns = _load_rtm_notebook()
    i_ext_vec, _, f_vec_periodic = ns.compute_rtm_f_i_curves_2()
    assert len(f_vec_periodic) == len(i_ext_vec)
    assert f_vec_periodic[-1] > 0.0
