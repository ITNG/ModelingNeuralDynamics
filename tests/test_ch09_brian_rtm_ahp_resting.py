from pathlib import Path

from matlab_ref import load_notebook_as_module, run_matlab_script, trace_rmse

ROOT = Path(__file__).resolve().parents[1]


def test_rtm_ahp_resting_matches_matlab():
    ns = load_notebook_as_module(ROOT / "brian" / "chapter09.ipynb")
    import neurodynex3.tools.input_factory as input_factory

    current = input_factory.get_step_current(0, 600, ns.b2.ms, 0.0 * ns.b2.uA)
    sm = ns.simulate_RTM_AHP_neuron(current, 600 * ns.b2.ms)
    t = sm.t / ns.b2.ms
    ca = sm.ca[0]

    ref = run_matlab_script(ROOT / "matlab" / "09" / "RTM_AHP_RESTING", "make_figure.m", ["t", "ca"])

    rmse = trace_rmse(ref["t"], ref["ca"], t, ca)
    assert rmse < 0.001, f"ca RMSE vs MATLAB too high: {rmse:.6f}"
