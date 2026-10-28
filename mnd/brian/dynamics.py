"""Brian2 trajectories for phase-plane, bifurcation and population-rate chapters.

States are dimensionless numbers (voltage numbers denote mV); one source
model time unit is represented by a Brian millisecond. Normal forms and FHN
use abstract model time, so their numerical time axis is not a physiological
calibration. Every integration is performed by Brian2, never SciPy.
"""
from types import SimpleNamespace

import brian2 as b2
import numpy as np

from .models import MODEL_DEFAULTS, _RATES


def _simulate(equations, initial, parameters, duration, dt, *, threshold=None,
              reset=None, refractory=False):
    """Integrate an explicit isolated network and include both endpoints."""
    if duration <= 0 or dt <= 0:
        raise ValueError('duration and dt must be positive')
    size = max(np.atleast_1d(value).size for value in initial.values())
    group = b2.NeuronGroup(size, equations, method='rk4', dt=dt*b2.ms,
                           namespace=parameters, threshold=threshold,
                           reset=reset, refractory=refractory)
    for name, value in initial.items():
        setattr(group, name, value)
    monitor = b2.StateMonitor(group, list(initial), record=True)
    objects = [group, monitor]
    spikes = None
    if threshold is not None:
        spikes = b2.SpikeMonitor(group)
        objects.append(spikes)
    network = b2.Network(*objects)
    network.run(duration*b2.ms, namespace={})
    monitor.record_single_timestep()
    result = {'t': np.asarray(monitor.t/b2.ms)}
    result.update({name: np.asarray(getattr(monitor, name)) for name in initial})
    result['spike_t'] = np.array([]) if spikes is None else np.asarray(spikes.t/b2.ms)
    result['spike_i'] = np.array([],dtype=int) if spikes is None else np.asarray(spikes.i)
    return SimpleNamespace(**result)


def simulate_fhn(*, a=1.25, tau_n=15.625, current=-.5,
                 initial=(-1.,-1.), duration=400., dt=.01):
    """FitzHugh–Nagumo voltage/recovery trajectories, chapters 10 and 15."""
    return _simulate('''dv/dt=(v-v**3/3-n+drive)/ms : 1
                        dn/dt=(a*v-n)/(tau_n*ms) : 1''',
                     dict(v=initial[0],n=initial[1]),
                     dict(a=a,tau_n=tau_n,drive=current),duration,dt)


def simulate_saddle_node(*, a=.45, initial=(.7,.2),duration=200.,dt=.01):
    """Saturating planar saddle-node example from chapter 11 (b=1)."""
    return _simulate('''dx/dt=(-a*x+y)/ms : 1
                        dy/dt=(x**2/(1+x**2)-y)/ms : 1''',
                     dict(x=initial[0],y=initial[1]),dict(a=a),duration,dt)


def simulate_hopf(*, current=.5, kind='supercritical',initial=(.1,0.),
                  duration=30.,dt=.01):
    """Cartesian Hopf flow; subcritical case includes quintic saturation.

    Angular frequency is 1 radian per abstract source time unit, matching the
    chapter's spiral helpers. The radial laws are I*r-r^3 and I*r+r^3-r^5.
    """
    if kind not in ('supercritical','subcritical'):
        raise ValueError('kind must be supercritical or subcritical')
    radial='drive-r2' if kind=='supercritical' else 'drive+r2-r2**2'
    return _simulate(f'''dx/dt=(({radial})*x-y)/ms : 1
                         dy/dt=(({radial})*y+x)/ms : 1
                         r2=x**2+y**2 : 1''',
                     dict(x=initial[0],y=initial[1]),dict(drive=current),duration,dt)


def simulate_reduced(model='hh',*,current=10.,initial=(-70.,.1),duration=50.,dt=.01):
    """Two-state HH (h=.83-n) or RTM (h=1-n), chapters 10,12,14.

    Rates are in inverse source-time units. exprel evaluates removable
    singularities without introducing special cases in the voltage sweep.
    """
    if model not in ('hh','rtm'):
        raise ValueError('reduced model must be hh or rtm')
    parameters={**MODEL_DEFAULTS[model],
                'h_offset': .83 if model=='hh' else 1.,'drive':current}
    rates=_RATES[model]
    equations='''dv/dt=(g_na*m_inf**3*(h_offset-n)*(e_na-v)
                        +g_k*n**4*(e_k-v)+g_l*(e_l-v)+drive)/ms : 1
                 dn/dt=(alpha_n*(1-n)-beta_n*n)/ms : 1
                 m_inf=alpha_m/(alpha_m+beta_m) : 1
                 '''+rates
    return _simulate(equations,dict(v=initial[0],n=initial[1]),parameters,duration,dt)


def simulate_self_exciting_theta(*,current=-.05,w_max=.2,tau_w=20.,
                                 duration=50.,dt=.01):
    """Chapter 16 feedback: phase crossing sets w=w_max (not an increment)."""
    return _simulate('''dtheta/dt=(1-cos(theta)+(drive+w)*(1+cos(theta)))/ms : 1
                        dw/dt=-w/(tau_w*ms) : 1''',
                     dict(theta=np.pi/2,w=0.),dict(drive=current,w_max=w_max,tau_w=tau_w),
                     duration,dt,threshold='theta>pi',reset='theta=-pi;w=w_max')


def simulate_wilson_cowan(*,w_ee=1.5,w_ie=1.,w_ei=1.,w_ii=0.,
                          tau_e=5.,tau_i=10.,drive_e=20.,drive_i=0.,
                          initial=(50.,10.),duration=300.,dt=.01):
    """Chapter 22 rectified, saturating E–I population-rate ODEs."""
    equations='''dE/dt=(100*u_e**2/(900+u_e**2)-E)/(tau_e*ms) : 1
                 dI/dt=(100*u_i**2/(400+u_i**2)-I)/(tau_i*ms) : 1
                 u_e=clip(w_ee*E-w_ie*I+drive_e,0,inf) : 1
                 u_i=clip(w_ei*E-w_ii*I+drive_i,0,inf) : 1'''
    return _simulate(equations,dict(E=initial[0],I=initial[1]),
                     dict(w_ee=w_ee,w_ie=w_ie,w_ei=w_ei,w_ii=w_ii,
                          tau_e=tau_e,tau_i=tau_i,drive_e=drive_e,drive_i=drive_i),duration,dt)
