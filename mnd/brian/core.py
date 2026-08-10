import brian2 as b2
import numpy as np


def get_step_current(t_start, t_end, unit_time, amplitude, append_zero=True):
    """Create a step current as a Brian2 TimedArray.

    Port of neurodynex3.tools.input_factory.get_step_current (the
    EPFL/BMNN course's Brian2 helper package), reimplemented here so this
    repo doesn't depend on it for one function.

    Args:
        t_start (int): start of the step
        t_end (int): end of the step
        unit_time (Brian2 unit): unit of t_start and t_end, e.g. brian2.ms
        amplitude (Quantity): amplitude of the step, e.g. 3.5*brian2.uamp
        append_zero (bool, optional): if True, 0 amp is appended at
            t_end+1. Without that trailing 0, Brian2 reads out the last
            value in the array (=amplitude) for all indices > t_end.

    Returns:
        TimedArray: Brian2 TimedArray
    """
    assert isinstance(t_start, int), "t_start must be of type int"
    assert isinstance(t_end, int), "t_end must be of type int"
    assert b2.units.fundamentalunits.have_same_dimensions(amplitude, b2.amp), \
        "amplitude must have the dimension of current, e.g. brian2.uamp"

    size = t_end + 2 if append_zero else t_end + 1
    values = np.zeros((size, 1)) * b2.amp
    values[t_start: t_end + 1, 0] = amplitude
    return b2.TimedArray(values, dt=1.0 * unit_time)
