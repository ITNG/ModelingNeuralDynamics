# Progress: MATLAB -> Python/Brian2 port status

Unit of tracking is the MATLAB sub-example (one leaf folder = one figure/script).
Process order: chapter number ascending, then sub-example name within a chapter.

- [x] done and verified (visual + numeric vs MATLAB)
**Totals: python 85/256 (33%). brian: not yet audited at sub-example grain, tracked here going forward.**

- [ ] not done

## Chapter 01 - Modeling a Single Neuron

| sub-example | python | brian |
|---|---|---|
| HH_VOLTAGE_TRACE | [x] | [ ] |

## Chapter 03 - The Classical HH ODEs

| sub-example | python | brian |
|---|---|---|
| HH_GATING_VARIABLES | [x] | [ ] |

## Chapter 04 - Numerical Solution of HH ODEs

| sub-example | python | brian |
|---|---|---|
| HH_LIMIT_CYCLE | [x] | [ ] |
| HH_REFRACTORINESS | [x] | [ ] |
| HH_SOLUTION | [x] | [ ] |

## Chapter 05 - The Simple Model of Neurons in Rodent Brains

| sub-example | python | brian |
|---|---|---|
| ERISIR_VOLTAGE_TRACE | [x] | [ ] |
| ERISIR_VOLTAGE_TRACE_2 | [x] | [ ] |
| RTM_VOLTAGE_TRACE | [x] | [ ] |
| THREE_MODELS_GATING_VARIABLES | [x] | [ ] |
| WB_VOLTAGE_TRACE | [x] | [ ] |

## Chapter 07 - Linear Integrate and Fire (LIF) Neurons

| sub-example | python | brian |
|---|---|---|
| LIF_NEURON_WITH_HH | [x] | [ ] |
| LIF_VOLTAGE_TRACE | [x] | [ ] |
| LIF_VOLTAGE_TRACE_2 | [x] | [ ] |
| SUBTHR_FOR_HH | [x] | [ ] |
| TAU_M_FOR_HH | [x] | [ ] |

## Chapter 08 - Quadratic Integrate and Fire (QIF) and Theta Neurons

| sub-example | python | brian |
|---|---|---|
| QIF_INFINITE_THRESHOLD | [x] | [ ] |
| QIF_VOLTAGE_TRACE | [x] | [ ] |
| THETA_FIRING | [x] | [ ] |
| THREE_CIRCLES | [x] | [ ] |

## Chapter 09 - Spike Frequency Adaptation

| sub-example | python | brian |
|---|---|---|
| ADAPTATION_MAP | [x] | [ ] |
| CALCIUM_RISE | [x] | [ ] |
| LIF_ADAPT | [x] | [ ] |
| M_CURRENT | [x] | [ ] |
| RTM_AHP | [x] | [ ] |
| RTM_AHP_RESTING | [x] | [ ] |
| RTM_M | [x] | [ ] |
| RTM_M_RESTING | [x] | [ ] |
| V_V_TILDE | [x] | [ ] |

## Chapter 10 - The Slow Fast Phase Plane

| sub-example | python | brian |
|---|---|---|
| FN | [x] | [ ] |
| HH_CYCLE_SPEED | [x] | [ ] |
| HH_H_PLUS_N | [x] | [ ] |
| HH_NULLCLINES_PLUS_SOLUTION | [x] | [ ] |
| REDUCED_HH | [x] | [ ] |

## Chapter 11 - The Saddle Node Bifurcation

| sub-example | python | brian |
|---|---|---|
| SADDLE_NODE_BIFURCATION | [x] | [ ] |

## Chapter 12 - Two Dimensional Bifurcation Analysis

| sub-example | python | brian |
|---|---|---|
| RTM_2D_FP | [x] | [ ] |
| RTM_2D_INVARIANT_CYCLE | [x] | [ ] |

## Chapter 13 - Hopf Bifurcations

| sub-example | python | brian |
|---|---|---|
| HOPF_SUB | [x] | [ ] |
| HOPF_SUB_2 | [x] | [ ] |
| HOPF_SUB_BIF_DIAG | [x] | [ ] |
| HOPF_SUB_BIF_DIAG_2 | [x] | [ ] |
| HOPF_SUB_PHASE_PLANE | [x] | [ ] |
| HOPF_SUB_PHASE_PLANE_2 | [x] | [ ] |
| HOPF_SUP | [x] | [ ] |
| HOPF_SUP_BIF_DIAG | [x] | [ ] |
| HOPF_SUP_PHASE_PLANE | [x] | [ ] |

## Chapter 14

| sub-example | python | brian |
|---|---|---|
| ERISIR_2D_FP | [ ] | [ ] |
| ERISIR_REDUCED | [ ] | [ ] |
| HH_REDUCED_COUNT_FP | [ ] | [ ] |
| HH_REDUCED_CYCLE_DISTANCE | [ ] | [ ] |
| HH_REDUCED_FIXED_POINTS | [ ] | [ ] |
| HH_REDUCED_FP_EVS | [ ] | [ ] |
| HH_REDUCED_REPELLING_CYCLE | [ ] | [ ] |

## Chapter 15

| sub-example | python | brian |
|---|---|---|
| CANARD | [ ] | [ ] |
| CANARD_2 | [ ] | [ ] |
| FITZHUGH_NAGUMO_MACRO | [ ] | [ ] |
| FITZHUGH_NAGUMO_MICRO | [ ] | [ ] |
| HH_REDUCED_BIF_DIAG | [ ] | [ ] |
| MMOS | [ ] | [ ] |

## Chapter 16

| sub-example | python | brian |
|---|---|---|
| INAPIK_FIXED_POINTS | [ ] | [ ] |
| INAPIK_PHASE_PLANE | [ ] | [ ] |
| SELF_EXCITING_THETA_NEURON | [ ] | [ ] |
| SELF_EXCITING_THETA_SMOOTH | [ ] | [ ] |
| SETN_PHASE_PLANE | [ ] | [ ] |

## Chapter 17 - Frequency Current Curves

| sub-example | python | brian |
|---|---|---|
| ERISIR_F_I_CURVE | [ ] | [ ] |
| HH_F_I_CURVE | [x] | [ ] |
| HH_REDUCED_F_I_CURVE | [ ] | [ ] |
| INAPIK_F_I_CURVE | [ ] | [ ] |
| INAPIK_SADDLE_CYCLE_DISTANCE | [ ] | [ ] |
| LIF_F_I_CURVE | [ ] | [ ] |
| RTM_F_I_CURVE | [x] | [ ] |
| RTM_F_I_CURVE_AT_ONSET | [x] | [ ] |
| RTM_WITH_M_CURRENT_F_I | [ ] | [ ] |
| SETN_F_I | [ ] | [ ] |
| THETA_F_I_CURVE | [ ] | [ ] |
| WB_F_I_CURVE | [ ] | [ ] |
| WB_F_I_CURVE_AT_ONSET | [ ] | [ ] |

## Chapter 18

| sub-example | python | brian |
|---|---|---|
| ERISIR_BISTABLE | [ ] | [ ] |
| ERISIR_BISTABLE_GATES | [ ] | [ ] |
| ERISIR_BISTABLE_LIMITED_H | [ ] | [ ] |
| HH_BISTABLE | [ ] | [ ] |
| HH_BISTABLE_GATES | [ ] | [ ] |
| HH_BISTABLE_LIMITED_N | [ ] | [ ] |
| H_CURRENT | [ ] | [ ] |
| PLOT_MODIFIED_TAU_R | [ ] | [ ] |
| RTM_F_I_CURVE_WITH_I_H | [ ] | [ ] |
| RTM_VOLTAGE_TRACE_WITH_I_H | [ ] | [ ] |
| RTM_WITH_I_H_BISTABLE | [ ] | [ ] |
| RTM_WITH_I_H_BISTABLE_GATES | [ ] | [ ] |
| RTM_WITH_I_H_LIMITED_R | [ ] | [ ] |
| RTM_WITH_I_M_BISTABLE | [ ] | [ ] |
| RTM_WITH_I_M_BISTABLE_GATES | [ ] | [ ] |
| RTM_WITH_I_M_LIMITED_W | [ ] | [ ] |

## Chapter 19

| sub-example | python | brian |
|---|---|---|
| ELLIPSES | [ ] | [ ] |
| ERISIR_PLUS_SLOW_I_K | [ ] | [ ] |
| ERISIR_SHOW_SLOW_I_K | [ ] | [ ] |
| INAPIK_PLUS_SLOW_I_K | [ ] | [ ] |
| INAPIK_PLUS_SLOW_I_K_3D | [ ] | [ ] |
| INAPIK_PLUS_STRONG_SLOW_I_K | [ ] | [ ] |
| INAPIK_PLUS_WEAK_SLOW_I_K | [ ] | [ ] |
| INAPIK_SHOW_SLOW_I_K | [ ] | [ ] |
| SQUARE_WAVES | [ ] | [ ] |

## Chapter 20 - Chemical Synapses

| sub-example | python | brian |
|---|---|---|
| B_JAHR_STEVENS | [ ] | [ ] |
| RTM_PLOT_Q | [x] | [ ] |
| RTM_PLOT_S | [x] | [ ] |
| RTM_PLOT_S_PRESCRIBE_TAU_PEAK | [x] | [ ] |
| RTM_PLOT_S_TWO_VARIABLES | [x] | [ ] |
| RTM_WITH_AUTAPSE_F_I_CURVE | [ ] | [ ] |
| S_BUILDUP | [x] | [ ] |
| S_SLOW_BUILDUP | [x] | [ ] |

## Chapter 21 - Gap Junctions

| sub-example | python | brian |
|---|---|---|
| LIF_NETWORK_WITH_GJ | [ ] | [ ] |
| RESET_THRESHOLD | [ ] | [ ] |
| WB_NETWORK_WITH_GJ | [x] | [ ] |
| WB_NETWORK_WITH_GJ_SUBTHRESHOLD | [ ] | [ ] |

## Chapter 22 - A Wilson Cowan Model of an Oscillatory E-I Network

| sub-example | python | brian |
|---|---|---|
| WILSON_COWAN_E_AND_I | [x] | [ ] |
| WILSON_COWAN_LOWERING_W_EE | [ ] | [ ] |
| WILSON_COWAN_PHASE_PLANE | [x] | [ ] |
| WILSON_COWAN_RASTERGRAM | [x] | [ ] |

## Chapter 23

| sub-example | python | brian |
|---|---|---|
| LIF_ENTRAINMENT | [ ] | [ ] |
| PLOT_F_ENTRAINMENT | [ ] | [ ] |
| PLOT_F_ENTRAINMENT_2 | [ ] | [ ] |
| WB_ENTRAINMENT_INTERVALS | [ ] | [ ] |
| WB_NEURON_ENTRAINED | [ ] | [ ] |
| WB_NEURON_IRREGULAR | [ ] | [ ] |
| WB_NEURON_N_TO_ONE | [ ] | [ ] |

## Chapter 24 - Synchronization by Fast Recurrent Excitation

| sub-example | python | brian |
|---|---|---|
| RTM_E_TO_E_HETEROGENEOUS | [ ] | [ ] |
| RTM_E_TO_E_NETWORK_1 | [x] | [ ] |
| RTM_E_TO_E_NETWORK_2 | [ ] | [ ] |
| RTM_SPLAY | [x] | [ ] |
| RTM_SYNC | [ ] | [ ] |
| RTM_TWO_CELL_NETWORK | [ ] | [ ] |

## Chapter 25 - Phase Response Curves (PRCs)

| sub-example | python | brian |
|---|---|---|
| MISC_PRC | [ ] | [ ] |
| PHASE_SHIFT | [ ] | [ ] |
| RTM_INTERACTION_FUNCTION | [ ] | [ ] |
| RTM_PRC | [ ] | [ ] |
| RTM_PRC_SHORT | [ ] | [ ] |
| RTM_PRC_SHORT_AND_WEAK | [ ] | [ ] |
| RTM_PRC_THREE_WEAK_ONES | [ ] | [ ] |
| RTM_PRC_WEAK | [ ] | [ ] |
| THETA_F | [ ] | [ ] |
| THETA_F_TILDE | [ ] | [ ] |
| THETA_PRC | [ ] | [ ] |
| THETA_PRC_SHORT_WEAK | [ ] | [ ] |
| WB_PRC_INHIBITORY_PULSE | [ ] | [ ] |

## Chapter 26

| sub-example | python | brian |
|---|---|---|
| ABSTRACT_PULSE_COUPLING_1 | [ ] | [ ] |
| ABSTRACT_PULSE_COUPLING_2 | [ ] | [ ] |
| ABSTRACT_PULSE_COUPLING_3 | [ ] | [ ] |
| ABSTRACT_PULSE_COUPLING_4 | [ ] | [ ] |
| ABSTRACT_PULSE_COUPLING_5 | [ ] | [ ] |
| F_TILDE | [ ] | [ ] |
| RTM_PLOT_G | [ ] | [ ] |
| TWO_PULSE_COUPLED_OSC | [ ] | [ ] |
| TWO_PULSE_COUPLED_OSC_2 | [ ] | [ ] |

## Chapter 27

| sub-example | python | brian |
|---|---|---|
| THREE_DELAYED_PULSE_COUPLED_OSC | [ ] | [ ] |
| TWO_DELAYED_PULSE_COUPLED_OSC | [ ] | [ ] |
| TWO_THETA_NEURONS | [ ] | [ ] |

## Chapter 28

| sub-example | python | brian |
|---|---|---|
| PLOT_D_TWO_FIXED_POINTS | [ ] | [ ] |
| WEAKLY_COUPLED_1 | [ ] | [ ] |
| WEAKLY_COUPLED_2 | [ ] | [ ] |
| WEAKLY_COUPLED_HETEROGENEOUS_1 | [ ] | [ ] |

## Chapter 29

| sub-example | python | brian |
|---|---|---|
| ILLUSTRATE_P0_AND_P1 | [ ] | [ ] |
| LIF_CONDITION_NUMBERS | [ ] | [ ] |
| LIF_P_AND_S | [ ] | [ ] |
| LIF_WITH_INHIBITORY_PULSE | [ ] | [ ] |
| RIVER | [ ] | [ ] |
| RTM_CONDITION_NUMBERS | [ ] | [ ] |
| RTM_WITH_INHIBITORY_PULSE | [ ] | [ ] |

## Chapter 30 - The PING Model of Gamma Rhythms

| sub-example | python | brian |
|---|---|---|
| 2_CELL_PING | [x] | [ ] |
| 2_CELL_PING_CONDITION_NUMBERS | [x] | [ ] |
| PING_1 | [x] | [ ] |
| PING_2 | [x] | [ ] |
| PING_3 | [x] | [ ] |
| PING_4 | [x] | [ ] |
| PING_5 | [ ] | [ ] |
| PING_6 | [ ] | [ ] |
| PING_7 | [ ] | [ ] |
| PING_8 | [ ] | [ ] |
| PING_9 | [ ] | [ ] |

## Chapter 31 - ING Rhythms

| sub-example | python | brian |
|---|---|---|
| 1_CELL_ING | [x] | [ ] |
| 1_CELL_ING_CONDITION_NUMBERS | [ ] | [ ] |
| ABSTRACT_PULSE_COUPLING_INH | [ ] | [ ] |
| ABSTRACT_PULSE_COUPLING_INH_2 | [ ] | [ ] |
| ING_1 | [x] | [ ] |
| ING_10 | [ ] | [ ] |
| ING_2 | [ ] | [ ] |
| ING_3 | [ ] | [ ] |
| ING_4 | [ ] | [ ] |
| ING_5 | [ ] | [ ] |
| ING_6 | [ ] | [ ] |
| ING_7 | [ ] | [ ] |
| ING_8 | [ ] | [ ] |
| ING_9 | [ ] | [ ] |
| ING_ENTRAINING_E_CELLS | [ ] | [ ] |
| ING_ENTRAINING_E_CELLS_2 | [ ] | [ ] |

## Chapter 32

| sub-example | python | brian |
|---|---|---|
| M_CURRENT_PING_1 | [ ] | [ ] |
| M_CURRENT_PING_1_CLOSEUP | [ ] | [ ] |
| M_CURRENT_PING_1_FROM_REST | [ ] | [ ] |
| M_CURRENT_PING_2_CLOSEUP | [ ] | [ ] |
| M_CURRENT_PING_3_CLOSEUP | [ ] | [ ] |
| PING_CLUSTERS | [ ] | [ ] |
| PLOT_PHI | [ ] | [ ] |
| PLOT_PSI | [ ] | [ ] |
| PLOT_PSI_PHI | [ ] | [ ] |
| POISSON_PING_1 | [ ] | [ ] |
| POISSON_PING_2 | [ ] | [ ] |
| POISSON_PING_3 | [ ] | [ ] |
| POISSON_PING_3_VOLTAGE_TRACE | [ ] | [ ] |

## Chapter 33

| sub-example | python | brian |
|---|---|---|
| M_CURRENT_BETA_WITH_GJ | [ ] | [ ] |
| M_CURRENT_PING_4 | [ ] | [ ] |
| M_CURRENT_PING_5 | [ ] | [ ] |
| M_CURRENT_PING_6 | [ ] | [ ] |
| M_CURRENT_PING_7 | [ ] | [ ] |
| M_CURRENT_PING_8 | [ ] | [ ] |
| PINB_1 | [ ] | [ ] |
| PINB_2 | [ ] | [ ] |
| PINB_3 | [ ] | [ ] |

## Chapter 34 - Nested Gamma Theta Rhythms

| sub-example | python | brian |
|---|---|---|
| A_CURRENT | [x] | [ ] |
| EIO_1 | [x] | [ ] |
| OLM_WITH_H_AND_A_CURRENTS | [x] | [ ] |
| OLM_WITH_H_CURRENT | [x] | [ ] |
| PING_WITH_THETA_DRIVE | [x] | [ ] |
| PING_WITH_THETA_INHIBITION | [x] | [ ] |
| PRE_OLM_VOLTAGE_TRACE | [x] | [ ] |
| PRE_OLM_X_INF_TAU_X | [x] | [ ] |

## Chapter 35

| sub-example | python | brian |
|---|---|---|
| OSCILLATIONS | [ ] | [ ] |
| PERIODIC_INHIBITION | [ ] | [ ] |
| PERIODIC_INHIBITION_2 | [ ] | [ ] |
| PERIODIC_INHIBITION_3 | [ ] | [ ] |
| PERIODIC_INHIBITION_F_I_CURVE | [ ] | [ ] |
| PERIODIC_INHIBITION_F_I_CURVE_2 | [ ] | [ ] |
| RTM_F_I_CURVE_WITH_INHIBITION | [ ] | [ ] |
| RTM_F_I_CURVE_WITH_INHIBITION_2 | [ ] | [ ] |

## Chapter 36

| sub-example | python | brian |
|---|---|---|
| IDEALIZED_F_I_CURVE | [ ] | [ ] |
| RTM_F_I_CURVE_PULSED_EXCITATION | [ ] | [ ] |
| RTM_F_I_CURVE_PULSED_EXCITATION_2 | [ ] | [ ] |
| SQUARE_PULSES | [ ] | [ ] |

## Chapter 37

| sub-example | python | brian |
|---|---|---|
| NO_RESET | [ ] | [ ] |
| PING_THR_1 | [ ] | [ ] |
| PING_THR_1_ZOOM | [ ] | [ ] |
| THRESHOLDING | [ ] | [ ] |

## Chapter 38

| sub-example | python | brian |
|---|---|---|
| GAMMA_COHERENCE_1 | [ ] | [ ] |
| GAMMA_COHERENCE_2 | [ ] | [ ] |
| POISSON_PING_3_MISMATCHED_PULSES | [ ] | [ ] |
| POISSON_PING_3_PLUS_GREEN | [ ] | [ ] |
| POISSON_PING_3_PLUS_PULSES | [ ] | [ ] |

## Chapter 39 - Short-Term Depression and Facilitation

| sub-example | python | brian |
|---|---|---|
| PULSES | [x] | [ ] |
| RTM_WITH_DEPRESSING_AND_FACILITATING_S | [x] | [ ] |
| RTM_WITH_DEPRESSING_S | [x] | [ ] |
| WB_WITH_DEPRESSING_S | [ ] | [ ] |

## Chapter 40 - Spike Timing-Dependent Plasticity(STDP)

| sub-example | python | brian |
|---|---|---|
| ABBOTT_SONG | [x] | [ ] |
| PING_WITH_STDP | [ ] | [ ] |
| RTM_VOLTAGE_TRACE_WITH_A | [x] | [ ] |
| THREE_CELL_PING_1 | [x] | [ ] |
| THREE_CELL_PING_2 | [x] | [ ] |
| THREE_CELL_PING_3 | [x] | [ ] |
| THREE_CELL_PING_4 | [x] | [ ] |
| THREE_CELL_PING_5 | [ ] | [ ] |

