"""Numerical contracts for the Brian2 geometry and rate-model examples."""
import importlib
import numpy as np
import pytest
from scipy.integrate import solve_ivp

b2 = pytest.importorskip('brian2')
pytestmark = pytest.mark.brian2


@pytest.fixture(autouse=True)
def isolated_numpy_backend():
    previous = b2.prefs.codegen.target
    b2.prefs.codegen.target = 'numpy'
    try:
        yield
    finally:
        b2.prefs.codegen.target = previous


def dynamics():
    return importlib.import_module('mnd.brian.dynamics')


def reference(rhs, initial, t):
    return solve_ivp(rhs, (0, float(t[-1])), initial, t_eval=t,
                     rtol=1e-10, atol=1e-11).y


def test_fhn_recovery_timescale_and_drive_match_scipy():
    d = dynamics()
    out = d.simulate_fhn(duration=20, dt=.02, initial=(-1, -1))
    expected = reference(lambda t, y: [y[0]-y[0]**3/3-y[1]-.5,
                                      (1.25*y[0]-y[1])/15.625], [-1,-1], out.t)
    np.testing.assert_allclose(np.array([out.v[0], out.n[0]]), expected, atol=2e-6)


def test_saddle_node_trajectories_match_source_equations():
    d = dynamics()
    out = d.simulate_saddle_node(a=.45, initial=(.7,.2), duration=10, dt=.02)
    expected = reference(lambda t,y: [-.45*y[0]+y[1], y[0]**2/(1+y[0]**2)-y[1]],
                         [.7,.2],out.t)
    np.testing.assert_allclose(np.array([out.x[0],out.y[0]]),expected,atol=1e-7)


def test_hopf_supercritical_radius_approaches_predicted_cycle():
    d = dynamics()
    out=d.simulate_hopf(current=.5, initial=(.1,0), duration=30, dt=.02)
    radius=np.hypot(out.x[0],out.y[0])
    assert radius[-1] == pytest.approx(np.sqrt(.5),abs=2e-6)
    assert radius[0] == pytest.approx(.1)
    angle=np.unwrap(np.arctan2(out.y[0],out.x[0]))
    assert angle[-1] == pytest.approx(30.,abs=1e-5)


def test_quintic_hopf_two_basins_are_preserved():
    d=dynamics()
    out=d.simulate_hopf(current=-.2,kind='subcritical',initial=([.2,.8],[0,0]),
                        duration=40,dt=.02)
    radius=np.hypot(out.x[:,-1],out.y[:,-1])
    assert radius[0] < .001
    assert radius[1] == pytest.approx(np.sqrt(.5+np.sqrt(.05)),abs=2e-6)


@pytest.mark.parametrize('model,h_offset,gna,gk,gl,ena,ek,el',[
    ('hh',.83,120,36,.3,45,-82,-59),('rtm',1.,100,80,.1,50,-100,-67)])
def test_reduced_conductance_trajectories_match_independent_rates(model,h_offset,gna,gk,gl,ena,ek,el):
    d=dynamics()
    def rhs(t,y):
        v,n=y
        if model=='hh':
            am=(v+45)/10/(1-np.exp(-(v+45)/10));bm=4*np.exp(-(v+70)/18)
            an=.01*(-60-v)/(np.exp((-60-v)/10)-1);bn=.125*np.exp(-(v+70)/80)
        else:
            am=.32*(v+54)/(1-np.exp(-(v+54)/4));bm=.28*(v+27)/(np.exp((v+27)/5)-1)
            an=.032*(v+52)/(1-np.exp(-(v+52)/5));bn=.5*np.exp(-(v+57)/40)
        m=am/(am+bm)
        return [gna*m**3*(h_offset-n)*(ena-v)+gk*n**4*(ek-v)+gl*(el-v)+1,
                an*(1-n)-bn*n]
    out=d.simulate_reduced(model,current=1,initial=(-70,.05),duration=2,dt=.0005)
    np.testing.assert_allclose(np.array([out.v[0],out.n[0]]),
                               reference(rhs,[-70,.05],out.t),atol=2e-5)


def test_theta_self_excitation_sets_feedback_at_spikes():
    d=dynamics()
    out=d.simulate_self_exciting_theta(duration=50,dt=.01)
    assert len(out.spike_t)>0
    assert np.max(out.w) == pytest.approx(.2)
    assert np.min(out.theta)>=-np.pi
    assert np.max(out.theta)<=np.pi+.03


def test_wilson_cowan_rectified_saturating_rates_match_scipy():
    d=dynamics()
    def rhs(t,y):
        e,i=y;a=max(1.5*e-i+20,0);c=max(e,0)
        return [(100*a*a/(900+a*a)-e)/5,(100*c*c/(400+c*c)-i)/10]
    out=d.simulate_wilson_cowan(duration=40,dt=.02)
    np.testing.assert_allclose(np.array([out.E[0],out.I[0]]),
                               reference(rhs,[50,10],out.t),atol=2e-6)
