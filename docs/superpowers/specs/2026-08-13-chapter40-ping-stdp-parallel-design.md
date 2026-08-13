# Chapter 40 PING-with-STDP Parallel Acceleration Design

## Goal

Accelerate the default `simulate_ping_with_stdp()` cell without changing its equations, integration timestep, random initialization, public interface, or returned data.

## Design

The dense STDP derivative updates each postsynaptic row independently. Keep its trace preparation serial, parallelize the 200 independent row updates with Numba `prange`, and retain the sequential predictor/corrector timestep order. The two-cell STDP example retains a serial wrapper around the same inlined row calculation, avoiding parallel-launch overhead for two rows. Temporarily cap an excessive Numba worker count at eight during the full-network simulation, restoring the caller's previous setting afterward; smaller configured thread counts remain unchanged.

Compile the initializer loop as `_rtm_init_population` and keep `rtm_init_population` as an array-normalizing Python wrapper, preserving its acceptance of Python sequences. Its loop and results remain the same; use Numba-compatible explicit integer and boolean dtypes.

## Validation

Before implementation, add tests requiring parallel compilation and nopython compilation of initialization. After implementation, compare the new STDP derivative directly against a serial reference calculation, check deterministic short-simulation outputs, run all Chapter 40 tests, and benchmark the default simulation with compilation and warm runtime reported separately.
