# Chapter 40 Numba Acceleration Design

## Goal

Accelerate `simulate_three_cell_ping` and `sweep_three_cell_ping_ee` without changing their public interfaces, return values, solver, or scientific interpretation.

## Design

Keep SciPy's adaptive `odeint` integration and compile the repeatedly invoked numerical right-hand side with Numba. Mark the RTM/WB channel-rate helpers and `derivative_three_cell_ping` with `@njit`. Replace the four `np.matmul(matrix, vector)` expressions in that derivative with equivalent `matrix @ vector` expressions, which Numba supports for these arrays.

The Python orchestration functions, `SimpleNamespace` results, spike detection, and sweep loop remain unchanged. Compilation cost is paid on the first simulation call; subsequent simulations and sweep points reuse the compiled specialization.

## Validation

Add a regression test that executes a short simulation, verifies that the derivative acquired a Numba nopython specialization, and checks expected output shapes and finite values. Run the Chapter 40 tests, compare compiled and original trajectories/spike times over representative weak and strong couplings, and benchmark warm execution.
