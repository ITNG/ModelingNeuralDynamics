# Task 3 Execution Report: Remaining Non-Brian Suite

Date: 2026-08-11

Worktree: `/home/ziaee/git/02_ITNG_REPOs/ModelingNeuralDynamics/.worktrees/perf-python-notebook-tests`

Baseline HEAD: `4595fb0 test: prevent full Python notebook execution`

Follow-up plan commit: `fca954d docs: plan measured non-brian suite optimization`
Corrected follow-up plan commit: `937681f docs: correct non-brian performance plan`
Evidence-gate correction commit: `803905b docs: align dense scan evidence gates`
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

## Complete Ranked Call-Duration Profile

The monolithic command produced no duration table, so this replacement ranks every completed call-duration line obtained by bounded batches and exact-node runs. Reproduced calls appear once at their three-run median. Ties retain chapter/file collection order. Lower bounds cannot be numerically interleaved with completed values and are listed separately after the table.

| Rank | Call | Timing | Evidence |
|---:|---|---:|---|
| 1 | ch24 heterogeneous network | 107.86s | exact node, one run |
| 2 | ch30 PING 5 | 94.17s | bounded batch, one run |
| 3 | ch32 M-current PING 1 from rest | 78.30s | bounded batch, one run |
| 4 | ch32 M-current PING 1 | 74.18s | bounded batch, one run |
| 5 | ch32 M-current PING 3 close-up | 72.46s | bounded batch, one run |
| 6 | ch35 periodic-inhibition f-I curve | 64.08s | bounded batch, one run |
| 7 | ch31 entraining E cells 1 | 42.33s | bounded batch, one run |
| 8 | ch32 Poisson PING 2 | 38.28s | bounded batch, one run |
| 9 | ch38 population statistics | 36.71s | bounded batch, one run |
| 10 | ch32 Poisson PING 3 | 35.34s | bounded batch, one run |
| 11 | ch32 Poisson PING 3 voltage | 34.75s | bounded batch, one run |
| 12 | ch32 PING clusters | 32.20s | bounded batch, one run |
| 13 | ch37 PING threshold zoom | 31.77s | bounded batch, one run |
| 14 | ch39 three-cell PING 5 | 29.41s | bounded batch, one run |
| 15 | ch33 PINB 1 | 28.24s | bounded batch, one run |
| 16 | ch14 HH reduced cycle distance | 27.06s | three-run median |
| 17 | ch33 PINB 3 | 26.68s | bounded batch, one run |
| 18 | ch33 M-current PING 7 | 25.92s | bounded batch, one run |
| 19 | ch33 M-current PING 6 | 25.84s | bounded batch, one run |
| 20 | ch33 M-current beta with gap junction | 25.81s | bounded batch, one run |
| 21 | ch33 PINB 2 | 24.55s | bounded batch, one run |
| 22 | ch30 PING 7 | 17.73s | exact node, one run |
| 23 | ch30 PING 8 | 16.37s | exact node, one run |
| 24 | ch30 PING 9 | 14.99s | exact node, one run |
| 25 | ch31 ING 10 | 14.88s | bounded batch, one run |
| 26 | ch31 ING 9 | 14.20s | bounded batch, one run |
| 27 | ch31 ING 6 | 13.88s | bounded batch, one run |
| 28 | ch31 ING 5 | 13.87s | bounded batch, one run |
| 29 | ch31 ING 2 | 13.74s | bounded batch, one run |
| 30 | ch31 ING 8 | 13.64s | bounded batch, one run |
| 31 | ch31 ING 4 | 13.58s | bounded batch, one run |
| 32 | ch31 ING 7 | 13.39s | bounded batch, one run |
| 33 | ch31 ING 3 | 13.14s | bounded batch, one run |
| 34 | ch37 PING threshold | 12.99s | bounded batch, one run |
| 35 | ch33 M-current PING 5 | 11.62s | bounded batch, one run |
| 36 | ch33 M-current PING 4 | 11.29s | bounded batch, one run |
| 37 | ch24 E-to-E network 1 | 9.98s | exact node, one run |
| 38 | ch24 splay-state network | 8.71s | exact node, one run |
| 39 | ch38 reproducible and active | 7.41s | bounded batch, one run |
| 40 | ch38 phase response, mismatched pulses | 6.78s | bounded batch, one run |
| 41 | ch18 h-current | 6.58s | bounded batch, one run |
| 42 | ch38 phase response, plus pulses | 6.58s | bounded batch, one run |
| 43 | ch18 modified tau-r | 6.55s | bounded batch, one run |
| 44 | ch23 f-entrainment plot | 6.38s | bounded batch, one run |
| 45 | ch09 notebook smoke | 6.24s | bounded batch, one run |
| 46 | ch17 HH legacy f-I | 5.35s | three-run median |
| 47 | ch01 HH voltage trace | 4.50s | bounded batch, one run |
| 48 | ch38 gamma coherence 2 | 4.38s | bounded batch, one run |
| 49 | ch20 RTM q plot | 4.34s | bounded batch, one run |
| 50 | ch31 one-cell ING condition numbers | 3.66s | bounded batch, one run |
| 51 | ch38 gamma coherence 1 | 2.93s | bounded batch, one run |
| 52 | ch12 RTM 2D fixed point | 2.71s | bounded batch, one run |
| 53 | ch39 WB depressing-s | 2.40s | bounded batch, one run |
| 54 | ch32 phi plot | 2.35s | bounded batch, one run |
| 55 | ch27 delayed pulse-coupled oscillators | 2.00s | bounded batch, one run |
| 56 | ch35 oscillations | 1.89s | bounded batch, one run |
| 57 | ch14 Erisir 2D fixed point | 1.88s | bounded batch, one run |
| 58 | ch37 no-reset | 1.88s | bounded batch, one run |
| 59 | ch36 idealized f-I | 1.87s | bounded batch, one run |
| 60 | ch17 RTM onset legacy f-I | 1.11s | three-run median |
| 61 | ch17 RTM legacy f-I | 0.79s | three-run median |
| 62 | ch22 phase plane | 0.72s | bounded batch, one run |
| 63 | ch32 psi-phi plot | 0.43s | bounded batch, one run |
| 64 | ch22 rastergram | 0.42s | bounded batch, one run |
| 65 | ch35 periodic inhibition 2 | 0.39s | bounded batch, one run |
| 66 | ch20 s buildup | 0.38s | bounded batch, one run |
| 67 | ch20 slow s buildup | 0.37s | bounded batch, one run |
| 68 | ch08 theta firing | 0.34s | bounded batch, one run |
| 69 | ch32 psi plot | 0.23s | bounded batch, one run |
| 70 | ch20 prescribed-tau plot | 0.13s | bounded batch, one run |
| 71 | loader repeated-load isolation | 0.13s | bounded batch, one run |
| 72 | ch20 two-stage s simulation | 0.12s | bounded batch, one run |
| 73 | ch35 periodic inhibition | 0.09s | bounded batch, one run |
| 74 | ch35 periodic inhibition 3 | 0.07s | bounded batch, one run |
| 75 | ch10 HH h+n | 0.07s | bounded batch, one run |
| 76 | ch05 Erisir voltage trace | 0.06s | bounded batch, one run |
| 77 | ch07 HH subthreshold | 0.05s | bounded batch, one run |
| 78 | ch10 reduced HH | 0.04s | bounded batch, one run |
| 79 | ch05 WB 1996 voltage trace | 0.03s | bounded batch, one run |
| 80 | ch08 QIF voltage trace | 0.03s | bounded batch, one run |
| 81 | ch05 WB voltage trace | 0.03s | bounded batch, one run |
| 82 | ch13 notebook smoke | 0.03s | bounded batch, one run |
| 83 | ch07 LIF voltage trace | 0.02s | bounded batch, one run |
| 84 | documentation relative links | 0.02s | bounded batch, one run |
| 85 | ch03 HH gating variables | 0.01s | bounded batch, one run |
| 86 | ch07 LIF voltage trace 2 | 0.01s | bounded batch, one run |
| 87 | ch22 Wilson-Cowan E+I | 0.01s | bounded batch, one run |
| 88 | ch38 shared-input case | 0.01s | bounded batch, one run |
| 89 | whole-notebook-loader regression guard | 0.01s | bounded batch, one run |

Cap/lower-bound ordering, all above the 107.86s completed maximum: PING 6 >115s; ING entraining E cells 2 >115s; both chapter-35 RTM inhibition calls >115s; chapter-36 pulsed excitation variant 2 >115s; chapter-36 pulsed excitation variant 1 >113s contextual lower bound. Additional context-only cap-active calls without attributable numerical lower bounds were chapter-32 M-current PING 1 close-up, M-current PING 2 close-up, and Poisson PING 1; chapter-33 M-current PING 8; chapter-35 periodic-inhibition f-I curve 2; chapter-37 thresholding; the final chapter-38 phase-statistics case; and chapter-40 PING with STDP.

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

Commits: `fca954d` (original plan), `937681f` (review-corrected plan), `803905b` (evidence-gate correction)

The plan contains seven implementation tasks with exact files, functions, tests, measured baselines, numerical assertions, and after-times. Tasks 1-5 use only measured smaller workloads or callable smoke entry points. Task 6 first removes import-time curve execution, then accepts cached Numba only when independent Python/JIT equivalence and fresh-cache cold cost both pass. Task 7 re-verifies Brian deselection and the bounded non-Brian suite.

## Self-Review

- **Spec coverage:** Brian deselection, exact required profiling commands, three-run stability, bounded remaining-suite coverage, hotspot classification, a measured follow-up plan, and artifact-only commits are all covered.
- **No-placeholder review:** searched the plan for `TBD`, `TODO`, `implement later`, `fill in details`, generic error-handling language, “Similar to,” and ellipsis placeholders; no unresolved implementation placeholder remains.
- **Interface consistency:** test calls match planned signatures (`t_final_run`, `run_connectivity_panels`, `run_drive_panels`, `run_smoke`, `i_values`, and `i_ext_values`). Return tuple orders are stated at each producer and consumed in the matching test snippet.
- **Scope review:** no Brian file appears in the plan's modification list. No production/test optimization was implemented in Task 3.
- **Evidence review:** after-times are tied to measured cold-process baselines or measured smaller-workload trials, not source appearance.

## Concerns

1. The monolithic suite cannot currently provide a trustworthy total duration or final all-tests result within a practical profiling window; bounded file/node execution is required until the dominant loops are reduced.
2. `test_inventory_matches_all_python_chapter_directories` is a reproducible pre-existing failure caused by the converted chapter-18 notebook coexisting with its legacy directory. It is outside this performance task.
3. Several cap-active nodes have lower bounds rather than completed timings. They are recorded but excluded from three-run stability claims.
4. The chapter-35/36 Numba step must pass the plan's independent Python/JIT three-point `rtol=0, atol=1e-10` comparisons and fresh-cache cold-cost gate before the full-grid speedup is accepted; otherwise Numba is removed for that script.

## Fix Round 1: Review Findings

Date: 2026-08-11

Plan fix commit: `937681f docs: correct non-brian performance plan`

### Changes made

1. Replaced the incomplete hotspot table with a complete descending profile of all 89 completed call-duration lines. Stable chapter-14/17 reproductions occur once at their medians. The six attributable cap/lower-bound calls are separate because they cannot be numerically interleaved with completed values. The plan maps all 54 retained calls (46 completed calls at or above 5s, two required sub-5s chapter-17 reproductions, and six cap/lower-bound calls) to an exact task or an evidence-based deferral.
2. Revised Task 6 so loading any targeted chapter-35/36 script creates no curve arrays and runs no scan. Tests call `compute_f_i_curve` or `compute_f_i_curves`; the unchanged full-grid call and plotting occur only below the main guard.
3. Replaced rebound-global Numba aliases with separately named Python functions and JIT dispatchers. Equivalence goes through public `use_numba=False` and `use_numba=True` paths. Cold validation uses a fresh `mktemp` cache directory and warm validation reuses it. The plan removes Numba independently from any script whose cold end-to-end time is not a net gain.
4. Changed PING 6 to build and simulate panel 1, then build and simulate panel 2, then build and simulate panel 3, preserving module RNG consumption.
5. Narrowed Task 5 to the only measured variant, M-current PING 1. Its sole interface is `run_smoke(t_final_run=t_final)` returning `(t_e_spikes, i_e_spikes, t_i_spikes, i_i_spikes, v_plot)`; no from-rest interface remains.
6. Explicitly deferred all four unmeasured M-current variants rather than extrapolating the PING-1 trial.
7. Removed the optional threshold-correction step. Final verification forbids numerical-threshold changes outside the exact assertions specified in Tasks 1-6.

No production or test optimization was implemented.

### Covering validation commands and exact output

Ranking consistency and retained-hotspot coverage:

```bash
python - <<'PY'
from pathlib import Path
import re
plan = Path('docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md').read_text()
report = Path('.superpowers/sdd/2026-08-11-python-notebook-test-loader-performance/task-3-report.md').read_text()
section = report.split('## Complete Ranked Call-Duration Profile', 1)[1].split('## Classification Evidence', 1)[0]
rows = []
for line in section.splitlines():
    match = re.match(r'\| (\d+) \| (.*?) \| ([0-9.]+)s \|', line)
    if match:
        rows.append((int(match.group(1)), match.group(2), float(match.group(3))))
assert [rank for rank, _, _ in rows] == list(range(1, len(rows) + 1))
assert all(a[2] >= b[2] for a, b in zip(rows, rows[1:]))
retained = [(name, timing) for _, name, timing in rows if timing >= 5.0]
assert len(retained) == 46
for required in ['28.24s', '26.68s', '25.92s', '27.06s median', '5.35s median', '1.11s median', '0.79s median']:
    assert required in plan, required
plan_dispositions = re.findall(r'^\| (\d+) \| .*? \| .*? \| (Task \d|Defer:)', plan, flags=re.M)
assert [int(rank) for rank, _ in plan_dispositions] == list(range(1, 55))
assert all(disposition.startswith(('Task ', 'Defer:')) for _, disposition in plan_dispositions)
print(f'ranking: {len(rows)} completed calls, ranks contiguous and non-increasing')
print(f'coverage: {len(retained)} calls >=5s plus 2 required sub-5s chapter-17 calls; 6 cap/lower-bound calls')
print(f'dispositions: {len(plan_dispositions)}/54 retained calls mapped to a task or explicit deferral')
PY
```

```text
ranking: 89 completed calls, ranks contiguous and non-increasing
coverage: 46 calls >=5s plus 2 required sub-5s chapter-17 calls; 6 cap/lower-bound calls
dispositions: 54/54 retained calls mapped to a task or explicit deferral
```

Placeholder scan:

```bash
rg -n 'TBD|TODO|implement later|fill in details|Similar to|\.\.\.|# existing|first_kernel|threshold corrections|simulate_from_rest|rest_duration|driven_duration|py\.f_vec|return \[simulate' \
  docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md
```

Output: none; exit code 1, meaning none of the forbidden placeholder, removed interface, import-global, or old RNG-order patterns remained.

Markdown fence balance:

```bash
python - <<'PY'
from pathlib import Path
for name in [Path('docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md'), Path('.superpowers/sdd/2026-08-11-python-notebook-test-loader-performance/task-3-report.md')]:
    count = sum(line.startswith('```') for line in name.read_text().splitlines())
    print(f'{name}: fences={count}; balanced={count % 2 == 0}')
PY
```

```text
docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md: fences=56; balanced=True
.superpowers/sdd/2026-08-11-python-notebook-test-loader-performance/task-3-report.md: fences=80; balanced=True
```

Interface/name consistency:

```bash
rg -n 'run_smoke|simulate_from_rest|rest_duration|driven_duration|compute_f_i_curve|compute_f_i_curves|f_vec_tonic|f_vec_periodic|f_vec_constant|f_vec_pulsed' \
  docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md
```

Exact conclusions from the matching lines: `run_smoke` has one signature and one five-array return order; none of `simulate_from_rest`, `rest_duration`, or `driven_duration` occurs; Task 6 consistently uses `compute_f_i_curve(i_values, use_numba)` for LIF and `compute_f_i_curves(i_ext_values, use_numba)` for RTM; legacy result names appear only in absence checks, returned-value assignments, and the documented main path.

Diff checks before the plan commit:

```bash
git diff --check
git diff --stat
git status --short
```

```text
 .../task-3-report.md                               | 129 ++++++---
 ...026-08-11-non-brian-python-suite-performance.md | 315 +++++++++++----------
 2 files changed, 261 insertions(+), 183 deletions(-)
 M .superpowers/sdd/2026-08-11-python-notebook-test-loader-performance/task-3-report.md
 M docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md
```

`git diff --check` emitted no output and exited 0.

### Fix-round self-review

- Ranking ranks all completed call timings numerically and keeps imprecise lower bounds separate.
- Every retained call has a task or explicit deferral; in particular the former omitted 94.17s, 38.28s, 28.24s, 26.68s, and 25.92s calls are covered.
- No Task 6 target performs a scan during `load_python_port`; Python/JIT equivalence does not depend on rebound globals.
- PING-6 panel construction/simulation stays sequential, and Task 5 contains no unmeasured variant or contradictory interface.
- Task 7 contains no optional numerical-contract changes.

### Fix-round concerns

The complete table is a bounded-run substitute for the unavailable monolithic `--durations=50` result; timings come from the exact commands already transcribed above and do not establish a full-suite pass. Context-only cap-active calls remain explicitly unranked and deferred because the bounded run never produced attributable call durations for them.

## Fix Round 2: Dense-Scan Target and Import Consistency

Date: 2026-08-11

Plan fix commit: `803905b docs: align dense scan evidence gates`

### Changes made

1. Removed chapter-35 `PERIODIC_INHIBITION_F_I_CURVE_2` from every Task 6 file list, callable/equivalence instruction, cold/warm command, commit command, and Task 7 optimized-set command. It remains solely in the explicit context-only deferral set because no attributable duration or lower bound was measured.
2. Added a five-row Task 6 evidence/gate table. Every target now names its completed duration or lower bound and an exact fresh-cache gate.
3. Added exact `from numba import njit` and `from numba.extending import register_jitable` imports. All decorators use `@register_jitable`; dispatchers use bare `njit`; the plan forbids an undefined bare `numba` module name.

No production or test file was modified.

### Covering validation commands and exact output

Target/defer disjointness, evidence/gate coverage, and Numba-name consistency:

```bash
python - <<'PY'
from pathlib import Path
import re
p = Path('docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md')
text = p.read_text()
task6 = text.split('### Task 6:', 1)[1].split('### Task 7:', 1)[0]
pre_tasks = text.split('### Task 1:', 1)[0]
targets = set(re.findall(r'^\| `(test_[^`]+)` \| ([^|]+) \| ([^|]+) \|$', task6, re.M))
target_names = {name for name, _, _ in targets}
expected_targets = {
    'test_periodic_inhibition_f_i_curve_structure',
    'test_rtm_f_i_curve_with_inhibition_structure',
    'test_rtm_f_i_curve_with_inhibition_2_structure',
    'test_rtm_f_i_curve_pulsed_excitation_structure',
    'test_rtm_f_i_curve_pulsed_excitation_2_structure',
}
assert target_names == expected_targets
for name, evidence, gate in targets:
    assert re.search(r'(\d+\.\d+s completed|>\d+s (exact cap|contextual lower bound))', evidence.strip()), (name, evidence)
    assert re.fullmatch(r'<30s and <\d+(?:\.\d+)?s', gate.strip()), (name, gate)
deferred = {
    'M_CURRENT_PING_1_CLOSEUP', 'M_CURRENT_PING_2_CLOSEUP', 'POISSON_PING_1',
    'M_CURRENT_PING_8', 'PERIODIC_INHIBITION_F_I_CURVE_2', 'THRESHOLDING',
    'FINAL_PHASE_STATISTICS', 'PING_WITH_STDP',
}
assert 'chapter-35 periodic-inhibition f-I curve 2' in pre_tasks
assert 'test_periodic_inhibition_f_i_curve_2_structure' not in task6
task6_scripts = {
    'PERIODIC_INHIBITION_F_I_CURVE', 'RTM_F_I_CURVE_WITH_INHIBITION',
    'RTM_F_I_CURVE_WITH_INHIBITION_2', 'RTM_F_I_CURVE_PULSED_EXCITATION',
    'RTM_F_I_CURVE_PULSED_EXCITATION_2',
}
assert not (task6_scripts & deferred)
assert task6.count('from numba import njit') == 1
assert task6.count('from numba.extending import register_jitable') == 1
assert '@numba.' not in task6
assert '@register_jitable' in task6
assert 'njit(cache=True)' in task6
print('target/defer overlap: none (5 Task 6 targets; 8 explicit context-only deferrals)')
print('measured gate coverage: 5/5 Task 6 targets')
print('Numba imports: PASS (njit and register_jitable each imported exactly once; no @numba usage)')
PY
```

```text
target/defer overlap: none (5 Task 6 targets; 8 explicit context-only deferrals)
measured gate coverage: 5/5 Task 6 targets
Numba imports: PASS (njit and register_jitable each imported exactly once; no @numba usage)
```

Placeholder scan:

```bash
rg -n 'TBD|TODO|implement later|fill in details|Similar to|\.\.\.|# existing|first_kernel|threshold corrections|simulate_from_rest|rest_duration|driven_duration|py\.f_vec|return \[simulate' \
  docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md
```

Output: none; exit code 1.

Markdown fence balance:

```bash
python - <<'PY'
from pathlib import Path
for p in [Path('docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md'), Path('.superpowers/sdd/2026-08-11-python-notebook-test-loader-performance/task-3-report.md')]:
    n = sum(line.startswith('```') for line in p.read_text().splitlines())
    print(f'{p}: fence lines={n}; balanced={n % 2 == 0}')
PY
```

```text
docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md: fence lines=58; balanced=True
.superpowers/sdd/2026-08-11-python-notebook-test-loader-performance/task-3-report.md: fence lines=94; balanced=True
```

Diff checks before the plan commit:

```bash
git diff --check
git diff --stat
git status --short
```

```text
 ...026-08-11-non-brian-python-suite-performance.md | 40 ++++++++++++++--------
 1 file changed, 25 insertions(+), 15 deletions(-)
 M docs/superpowers/plans/2026-08-11-non-brian-python-suite-performance.md
```

`git diff --check` emitted no output and exited 0.

### Fix-round self-review and concerns

Task 6 contains exactly five targets and no deferred script. Every target has an evidence row, a cold threshold, an independent Python/JIT equivalence path, and identifiers supplied by exact imports. The only standing concern is unchanged: `PERIODIC_INHIBITION_F_I_CURVE_2` cannot receive a performance task until an attributable measurement exists.
