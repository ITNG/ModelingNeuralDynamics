# Task 3 Execution Report: Remaining Non-Brian Suite

Date: 2026-08-11

Worktree: `/home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics/.worktrees/perf-python-notebook-tests`

Baseline HEAD: `4595fb0 test: prevent full Python notebook execution`

Follow-up plan commit: `fca954d docs: plan measured non-brian suite optimization`
Follow-up plan: `docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md`

## Outcome

- Explicit Brian collection selected zero tests: `95 deselected`.
- Default collection selected `227/450` tests and deselected `223`; Brian-marked tests remained excluded.
- The exact monolithic `pytest -q --durations=50` command did not finish within the available profiling window. It emitted progress dots only, remained alive after the controller interruption, and was observed at `51:12` elapsed before the exact pytest PID was stopped with `SIGINT`. This is a profiling constraint, not a test failure.
- Bounded chapter/file profiling found no simulation assertion failures. Commands were capped at 115 seconds and cap-active nodes were narrowed by exact node ID or source location.
- The final miscellaneous/docs batch produced `120 passed, 1 failed in 3.03s`. The one failure is a pre-existing documentation-inventory inconsistency: `python/chapter18.ipynb` makes chapter 18 “converted,” while the legacy `python/18_Bistability_Resulting_from_Rebound_Firing/` directory still exists. No fix was made because this task is measurement and planning only.
- No production or test file was modified.

## Exact Required Commands and Output

### Brian deselection

Command:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest --collect-only -q tests/test_*_brian_*.py
```

Output:

```text
no tests collected (95 deselected) in 5.84s
```

Pytest exit code was 5 because no test was selected. This matches the requirement.

### Default collection

Command:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest --collect-only -q
```

Final output:

```text
227/450 tests collected (223 deselected) in 3.87s
```

### Exact monolithic profile

Command:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q --durations=50
```

Output before interruption consisted only of pytest progress dots; pytest never emitted a duration table or pass/fail summary. A later process check showed:

```text
3922770       51:12 /home/ziaee/envs/mnd/bin/python -m pytest -q --durations=50
```

The exact process was stopped outside the sandbox because the abandoned sandbox wrapper could no longer signal it:

```bash
kill -INT 3922770
ps -p 3922770 -o pid=,etime=,args=
```

The `kill` command exited 0 and the final `ps` emitted no output (exit 1), confirming the process was gone.

### Chapter-14 three-run reproduction

Command, run three times:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch14_hh_reduced_cycle_distance.py::test_hh_reduced_cycle_distance_fixed_points_are_sane \
  --durations=1
```

Outputs:

```text
27.59s call tests/test_ch14_hh_reduced_cycle_distance.py::test_hh_reduced_cycle_distance_fixed_points_are_sane
1 passed in 27.85s

26.55s call tests/test_ch14_hh_reduced_cycle_distance.py::test_hh_reduced_cycle_distance_fixed_points_are_sane
1 passed in 26.81s

27.06s call tests/test_ch14_hh_reduced_cycle_distance.py::test_hh_reduced_cycle_distance_fixed_points_are_sane
1 passed in 27.33s
```

Median: 27.06s. Deviations from the median: +1.96%, -1.88%, and 0%. All are within 15%, so this hotspot is retained.

### Chapter-17 three-run reproduction

Command, run three times:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch17_legacy_f_i_curves.py \
  --durations=3
```

Outputs:

```text
Run 1
5.35s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_hh_f_i_curve_legacy
1.11s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_rtm_f_i_curve_at_onset_legacy
0.79s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_rtm_f_i_curve_legacy
3 passed in 7.52s

Run 2
5.20s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_hh_f_i_curve_legacy
1.09s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_rtm_f_i_curve_at_onset_legacy
0.78s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_rtm_f_i_curve_legacy
3 passed in 7.32s

Run 3
5.52s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_hh_f_i_curve_legacy
1.15s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_rtm_f_i_curve_at_onset_legacy
0.83s call tests/test_ch17_legacy_f_i_curves.py::test_chapter17_notebook_rtm_f_i_curve_legacy
3 passed in 7.76s
```

Medians and maximum deviations:

| Node | Median | Maximum deviation | Retained |
|---|---:|---:|---|
| HH legacy | 5.35s | 3.18% | yes |
| RTM onset legacy | 1.11s | 3.60% | yes |
| RTM legacy | 0.79s | 5.06% | yes |

## Bounded Remaining-Suite Commands and Output

Every command below used `timeout --signal=INT 115s`. Duration lines are copied from pytest output. A `124` result means the command reached the cap and pytest reported `KeyboardInterrupt`; it does not mean an assertion failed.

### Chapters 01-08

```bash
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch01_hh_voltage_trace.py tests/test_ch03_hh_gating_variables.py \
  tests/test_ch05_erisir_voltage_trace.py tests/test_ch05_rtm_gating_variables.py \
  tests/test_ch05_wb_voltage_trace.py tests/test_ch05_wb_voltage_trace_1996.py \
  tests/test_ch07_hh_subthreshold.py tests/test_ch07_lif_voltage_trace.py \
  tests/test_ch07_lif_voltage_trace_2.py tests/test_ch08_qif_voltage_trace.py \
  tests/test_ch08_theta_firing.py --durations=20
```

```text
4.50s call tests/test_ch01_hh_voltage_trace.py::test_chapter01_notebook_hh_voltage_trace
0.34s call tests/test_ch08_theta_firing.py::test_chapter08_notebook_theta_firing
11 passed in 5.85s
```

### Chapters 09-18 excluding the reproduced chapter-14/17 nodes

```bash
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_ch09_notebook_smoke.py tests/test_ch10_hh_h_plus_n.py \
  tests/test_ch10_reduced_hh.py tests/test_ch12_rtm_2d_fp.py \
  tests/test_ch13_notebook_smoke.py tests/test_ch14_erisir_2d_fp.py \
  tests/test_ch18_h_current.py tests/test_ch18_plot_modified_tau_r.py --durations=20
```

```text
6.58s call tests/test_ch18_h_current.py::test_h_current_matches_matlab
6.55s call tests/test_ch18_plot_modified_tau_r.py::test_plot_modified_tau_r_matches_matlab
6.24s call tests/test_ch09_notebook_smoke.py::test_chapter09_notebook_defines_working_examples
2.71s call tests/test_ch12_rtm_2d_fp.py::test_rtm_2d_fp_bifurcation_shape
1.88s call tests/test_ch14_erisir_2d_fp.py::test_erisir_2d_fp_bifurcation_shape
8 passed in 24.77s
```

### Chapters 20-24 and narrowing

The combined command passed the first 11 tests and capped while entering chapter 24:

```text
6.38s call tests/test_ch23_plot_f_entrainment.py::test_plot_f_entrainment_matches_matlab
4.34s call tests/test_ch20_rtm_plot_q.py::test_rtm_plot_q_smoke
11 passed, 1 warning in 115.03s
KeyboardInterrupt at mnd/core.py:47
```

Exact chapter-24 narrowing:

```bash
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch24_rtm_e_to_e_network_1.py --durations=3
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch24_rtm_splay.py --durations=3
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch24_rtm_e_to_e_heterogeneous.py --durations=3
```

```text
9.98s call tests/test_ch24_rtm_e_to_e_network_1.py::test_rtm_e_to_e_network_smoke
1 passed in 10.23s

8.71s call tests/test_ch24_rtm_splay.py::test_rtm_splay_smoke
1 passed in 8.97s

107.86s call tests/test_ch24_rtm_e_to_e_heterogeneous.py::test_rtm_e_to_e_heterogeneous_is_sane
1 passed in 108.11s
```

### Chapters 27-30 and narrowing

```text
94.17s call tests/test_ch30_ping_5.py::test_ping_5_structure
2.00s call tests/test_ch27_three_delayed_pulse_coupled_osc.py::test_three_delayed_pulse_coupled_osc_structure
2 passed in 115.03s
KeyboardInterrupt at python/30_The_PING_Model_of_Gamma_Rhythms/PING_6/main.py:37
```

Exact PING 6-9 outputs:

```text
PING 6: no tests ran in 115.01s; KeyboardInterrupt at PING_6/main.py:50
PING 7: 17.73s call; 1 passed in 17.99s
PING 8: 16.37s call; 1 passed in 16.63s
PING 9: 14.99s call; 1 passed in 15.25s
```

### Chapter 31 and narrowing

```text
14.20s test_ing_9_structure
13.88s test_ing_6_structure
13.87s test_ing_5_structure
13.74s test_ing_2_structure
13.64s test_ing_8_structure
13.58s test_ing_4_structure
13.39s test_ing_7_structure
13.14s test_ing_3_structure
3.66s test_1_cell_ing_condition_numbers_matches_matlab
9 passed in 115.02s
KeyboardInterrupt at ING_10/main.py:38
```

Continuation:

```text
42.33s call tests/test_ch31_ing_entraining_e_cells.py::test_ing_entraining_e_cells_structure
14.88s call tests/test_ch31_ing_2_to_10.py::test_ing_10_structure
2 passed in 115.02s
KeyboardInterrupt at ING_ENTRAINING_E_CELLS_2/main.py:38

Exact ING_ENTRAINING_E_CELLS_2: no tests ran in 115.03s; KeyboardInterrupt at main.py:32
```

### Chapter 32 and narrowing

First batch:

```text
74.18s call test_m_current_ping_1_structure
2.35s call test_plot_phi_structure
0.43s call test_plot_psi_phi_structure
0.23s call test_plot_psi_structure
4 passed in 114.99s
KeyboardInterrupt at M_CURRENT_PING_1_CLOSEUP/main.py:52
```

Continuations:

```text
78.30s call test_m_current_ping_1_from_rest_structure
1 passed in 115.02s
KeyboardInterrupt at M_CURRENT_PING_2_CLOSEUP/main.py:65

72.46s call test_m_current_ping_3_closeup_structure
32.20s call test_ping_clusters_structure
2 passed in 115.00s
KeyboardInterrupt at POISSON_PING_1/main.py:60

38.28s call test_poisson_ping_2_structure
35.34s call test_poisson_ping_3_structure
34.75s call test_poisson_ping_3_voltage_trace_structure
3 passed in 108.62s
```

### Chapter 33 and narrowing

```text
25.92s call test_m_current_ping_7_structure
25.84s call test_m_current_ping_6_structure
25.81s call test_m_current_beta_with_gj_structure
11.62s call test_m_current_ping_5_structure
11.29s call test_m_current_ping_4_structure
5 passed in 115.01s
KeyboardInterrupt at M_CURRENT_PING_8/main.py:82

28.24s call test_pinb_1_structure
26.68s call test_pinb_3_structure
24.55s call test_pinb_2_structure
3 passed in 79.74s
```

### Chapters 35-37

```text
64.08s call tests/test_ch35_periodic_inhibition.py::test_periodic_inhibition_f_i_curve_structure
1.89s call test_oscillations_structure
5 passed in 114.99s
KeyboardInterrupt at PERIODIC_INHIBITION_F_I_CURVE_2/main.py:19

Exact RTM f-I with inhibition: no tests ran in 115.01s; KeyboardInterrupt at main.py:30
Exact RTM f-I with inhibition 2: no tests ran in 114.99s; KeyboardInterrupt at main.py:30

Chapter 36 first batch: 1.87s idealized f-I; 2 passed in 115.00s;
KeyboardInterrupt at RTM_F_I_CURVE_PULSED_EXCITATION/main.py:33
Exact pulsed-excitation variant 2: no tests ran in 115.00s;
KeyboardInterrupt at RTM_F_I_CURVE_PULSED_EXCITATION_2/main.py:42

31.77s call tests/test_ch37_thresholding_in_ping.py::test_ping_thr_1_zoom_structure
12.99s call test_ping_thr_1_structure
1.88s call test_no_reset_structure
3 passed in 114.99s
KeyboardInterrupt at THRESHOLDING/main.py:40
```

### Chapters 38-40

```text
36.71s call tests/test_ch38_gamma_coherence.py::test_poisson_ping_population_statistics_match_matlab
7.41s call test_poisson_ping_is_reproducible_and_active
6.78s call test_poisson_ping_phase_response[POISSON_PING_3_MISMATCHED_PULSES-29.0]
6.58s call test_poisson_ping_phase_response[POISSON_PING_3_PLUS_PULSES-31.0]
4.38s call test_gamma_coherence_2_replaces_feedback_with_mean_inhibition
2.93s call test_gamma_coherence_1_produces_aligned_traces_and_spikes
8 passed, 2 deselected in 115.01s
KeyboardInterrupt at python/38_Gamma_Coherence/model.py:15

29.41s call tests/test_ch39_ch40_remaining.py::test_three_cell_ping_5_structure
2.40s call tests/test_ch39_ch40_remaining.py::test_wb_with_depressing_s_structure
2 passed in 114.98s
KeyboardInterrupt at PING_WITH_STDP/main.py:229
```

### Core, loader, and documentation

Command:

```bash
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_mnd_core.py tests/test_notebook_definition_loader.py \
  tests/test_python_chapter_docs.py --durations=30
```

Output:

```text
0.13s call tests/test_notebook_definition_loader.py::test_repeated_loads_reuse_compilation_but_isolate_namespaces
0.02s call tests/test_python_chapter_docs.py::test_all_documentation_relative_links_resolve
0.01s call tests/test_notebook_definition_loader.py::test_python_notebook_tests_never_use_whole_notebook_loader
FAILED tests/test_python_chapter_docs.py::test_inventory_matches_all_python_chapter_directories
1 failed, 120 passed in 3.03s
```

Exact reproduction:

```bash
/home/ziaee/envs/mnd/bin/python -m pytest -q \
  tests/test_python_chapter_docs.py::test_inventory_matches_all_python_chapter_directories
```

```text
Extra items in the right set:
'18_Bistability_Resulting_from_Rebound_Firing'
1 failed in 0.09s
```

### Exact bounded command log

The duration/summary output for each command is transcribed in the chapter sections above. These are the exact invocations, in execution order after chapters 01-18:

```bash
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch20_rtm_plot_q.py tests/test_ch20_rtm_plot_s_prescribe_tau_peak.py tests/test_ch20_rtm_plot_s_two_variables.py tests/test_ch20_s_buildup.py tests/test_ch20_s_slow_buildup.py tests/test_ch22_wilson_cowan_e_and_i.py tests/test_ch22_wilson_cowan_phase_plane.py tests/test_ch22_wilson_cowan_rastergram.py tests/test_ch23_plot_f_entrainment.py tests/test_ch23_plot_f_entrainment_2.py tests/test_ch24_rtm_e_to_e_heterogeneous.py tests/test_ch24_rtm_e_to_e_network_1.py tests/test_ch24_rtm_splay.py --durations=30

timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch24_rtm_e_to_e_network_1.py --durations=3
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch24_rtm_splay.py --durations=3
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch24_rtm_e_to_e_heterogeneous.py --durations=3

timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch27_three_delayed_pulse_coupled_osc.py tests/test_ch30_ping_5.py tests/test_ch30_ping_6_to_9.py --durations=20
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch30_ping_6_to_9.py::test_ping_6_structure --durations=3
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch30_ping_6_to_9.py::test_ping_7_structure --durations=3
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch30_ping_6_to_9.py::test_ping_8_structure --durations=3
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch30_ping_6_to_9.py::test_ping_9_structure --durations=3

timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch31_1_cell_ing_condition_numbers.py tests/test_ch31_ing_2_to_10.py tests/test_ch31_ing_entraining_e_cells.py --durations=30
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch31_ing_2_to_10.py::test_ing_10_structure tests/test_ch31_ing_entraining_e_cells.py --durations=10
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch31_ing_entraining_e_cells.py::test_ing_entraining_e_cells_2_structure --durations=3

timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch32_m_current_ping_poisson_ping.py --durations=30
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_1_from_rest_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_2_closeup_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_3_closeup_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_ping_clusters_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_1_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_2_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_3_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_3_voltage_trace_structure --durations=20
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch32_m_current_ping_poisson_ping.py::test_m_current_ping_3_closeup_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_ping_clusters_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_1_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_2_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_3_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_3_voltage_trace_structure --durations=20
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_2_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_3_structure tests/test_ch32_m_current_ping_poisson_ping.py::test_poisson_ping_3_voltage_trace_structure --durations=10

timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch33_m_current_ping_and_pinb.py --durations=30
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch33_m_current_ping_and_pinb.py::test_pinb_1_structure tests/test_ch33_m_current_ping_and_pinb.py::test_pinb_2_structure tests/test_ch33_m_current_ping_and_pinb.py::test_pinb_3_structure --durations=10

timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch35_periodic_inhibition.py --durations=20
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_structure tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_2_structure --durations=10
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch35_periodic_inhibition.py::test_rtm_f_i_curve_with_inhibition_2_structure --durations=3

timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch36_f_i_curves_pulsed_excitation.py --durations=15
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch36_f_i_curves_pulsed_excitation.py::test_rtm_f_i_curve_pulsed_excitation_2_structure --durations=3
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch37_thresholding_in_ping.py --durations=15
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch38_gamma_coherence.py --durations=25
timeout --signal=INT 115s /home/ziaee/envs/mnd/bin/python -m pytest -q tests/test_ch39_ch40_remaining.py --durations=10
```

## Ranked Confirmed Hotspots

Completed timings are ranked first. Stable three-run medians are marked `median`; single bounded-run measurements are marked `one run`. Exact cap results are lower bounds.

| Rank | Call | Timing | Evidence | Classification |
|---:|---|---:|---|---|
| 1 | ch30 `test_ping_6_structure` | >115s | exact cap | oversized test workload |
| 2 | ch31 `test_ing_entraining_e_cells_2_structure` | >115s | exact cap | oversized test workload |
| 3 | ch35 `test_rtm_f_i_curve_with_inhibition_structure` | >115s | exact cap | simulation loop |
| 4 | ch35 `test_rtm_f_i_curve_with_inhibition_2_structure` | >115s | exact cap | simulation loop |
| 5 | ch36 `test_rtm_f_i_curve_pulsed_excitation_2_structure` | >115s | exact cap | simulation loop |
| 6 | ch36 `test_rtm_f_i_curve_pulsed_excitation_structure` | >113s | narrowed batch lower bound | simulation loop |
| 7 | ch24 `test_rtm_e_to_e_heterogeneous_is_sane` | 107.86s | one run | oversized test workload |
| 8 | ch30 `test_ping_5_structure` | 94.17s | one run | simulation loop |
| 9 | ch32 `test_m_current_ping_1_from_rest_structure` | 78.30s | one run | oversized test workload |
| 10 | ch32 `test_m_current_ping_1_structure` | 74.18s | one run | oversized test workload |
| 11 | ch32 `test_m_current_ping_3_closeup_structure` | 72.46s | one run | oversized test workload |
| 12 | ch35 `test_periodic_inhibition_f_i_curve_structure` | 64.08s | one run | simulation loop |
| 13 | ch31 `test_ing_entraining_e_cells_structure` | 42.33s | one run | simulation loop |
| 14 | ch32 `test_poisson_ping_2_structure` | 38.28s | one run | simulation loop |
| 15 | ch38 `test_poisson_ping_population_statistics_match_matlab` | 36.71s | one run | oversized test workload |
| 16 | ch32 `test_poisson_ping_3_structure` | 35.34s | one run | simulation loop |
| 17 | ch32 `test_poisson_ping_3_voltage_trace_structure` | 34.75s | one run | simulation loop |
| 18 | ch32 `test_ping_clusters_structure` | 32.20s | one run | simulation loop |
| 19 | ch37 `test_ping_thr_1_zoom_structure` | 31.77s | one run | simulation loop |
| 20 | ch39/40 `test_three_cell_ping_5_structure` | 29.41s | one run | simulation loop |
| 21 | ch14 cycle-distance test | 27.06s | three-run median | oversized test workload |
| 22 | ch17 HH legacy f-I | 5.35s | three-run median | oversized test workload |
| 23 | ch17 RTM onset legacy f-I | 1.11s | three-run median | oversized test workload |
| 24 | ch17 RTM legacy f-I | 0.79s | three-run median | oversized test workload |

Additional cap-active calls with only contextual lower bounds were recorded but not ranked as completed timings: chapter-32 M-current PING 1 close-up, M-current PING 2 close-up, and Poisson PING 1; chapter-33 M-current PING 8; chapter-35 periodic-inhibition f-I curve 2; chapter-37 thresholding; the final chapter-38 phase-statistics case; and chapter-40 PING with STDP.

## Classification Evidence

- **Notebook loading:** no retained hotspot. The loader regression test took 0.13s, and all slow duration entries were pytest call time inside explicitly invoked simulations. Full-notebook execution remains prohibited by the passing regression guard.
- **Oversized test workload:** chapter 14 runs three currents, two directions, 100,000 steps per trajectory even though `t_final=200` preserves every assertion. Chapter 17 smoke tests use 200ms where 100ms preserves and strengthens the qualitative assertions. Chapter 24 uses 30 neurons for 200ms where a seeded 4-neuron/100ms trial preserves shape, ranges, zero diagonal, and repeated firing. PING 6 and ING entrainment use published figure horizons at import time; measured shortened horizons preserved their existing spike-count assertions.
- **Simulation loop:** the legacy chapter-30 through chapter-40 scripts spend their call time in midpoint/RK loops. The dense chapter-35/36 f-I tests have nested current-by-time loops and require their full 101/201-point grids for plateau assertions, making them the measured Numba candidates.
- **SciPy solver:** no SciPy solver reached the dominant retained set. The slower bifurcation calls (chapter 12 at 2.71s and chapter 14 Erisir at 1.88s) were below the planned threshold and were not selected for Numba.

## Measured Workload Trials Used by the Plan

These commands did not modify tracked files; they loaded definitions and called existing functions with explicit parameters.

| Function/workload | Output | Wall |
|---|---|---:|
| ch14 `simulate_hh_reduced_cycle_distance(t_final=200.0)` | all 3 attracting/repelling panels enter the window | 6.64s |
| ch17 three legacy functions at `t_final=100.0` | HH `[0,62.9167,75.1984]`; RTM `[0,0,42.6533]`; onset zeros | 3.99s total |
| ch24 heterogeneous `N=4,t_final=100,seed=63806` | counts `[2,2,2,2]`, shapes `(4,)`/`(4,4)` | 19.07s |
| ch30 PING 6 three panels, `t_final=200` | counts `(1888,437)`, `(1947,395)`, `(1943,424)` | 38.63s |
| ch31 ING E2 three drives, `t_final=100` | counts `(2059,510)`, `(2231,419)`, `(2341,447)` | 29.57s |
| ch32 M-current PING 1, `t_final=200` | 368 E, 282 I, max voltage 48.377mV | 21.81s |

## Follow-up Plan

Path: `docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md`

Commit: `fca954d`

The plan contains seven implementation tasks with exact files, functions, tests, measured baselines, numerical assertions, and after-times. Tasks 1-5 use smaller workloads or callable smoke entry points. Task 6 uses cached Numba only for dense chapter-35/36 scans whose grid resolution is part of the numerical contract. Task 7 re-verifies Brian deselection and the bounded non-Brian suite.

## Self-Review

- **Spec coverage:** Brian deselection, exact required profiling commands, three-run stability, bounded remaining-suite coverage, hotspot classification, a measured follow-up plan, and artifact-only commits are all covered.
- **No-placeholder review:** searched the plan for `TBD`, `TODO`, `implement later`, `fill in details`, generic error-handling language, “Similar to,” and ellipsis placeholders; no unresolved implementation placeholder remains.
- **Interface consistency:** test calls match planned signatures (`t_final_run`, `run_connectivity_panels`, `run_drive_panels`, `run_smoke`, `simulate_from_rest`, and `i_ext_values`). Return tuple orders are stated at each producer and consumed in the matching test snippet.
- **Scope review:** no Brian file appears in the plan's modification list. No production/test optimization was implemented in Task 3.
- **Evidence review:** after-times are tied to measured cold-process baselines or measured smaller-workload trials, not source appearance.

## Concerns

1. The monolithic suite cannot currently provide a trustworthy total duration or final all-tests result within a practical profiling window; bounded file/node execution is required until the dominant loops are reduced.
2. `test_inventory_matches_all_python_chapter_directories` is a reproducible pre-existing failure caused by the converted chapter-18 notebook coexisting with its legacy directory. It is outside this performance task.
3. Several cap-active nodes have lower bounds rather than completed timings. They are recorded but excluded from three-run stability claims.
4. The chapter-35/36 Numba step must pass the plan's three-point `rtol=0, atol=1e-10` comparisons before the full-grid speedup is accepted.
