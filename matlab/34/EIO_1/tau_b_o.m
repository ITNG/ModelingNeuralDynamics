function tau_b_o=tau_b_o(v)
tau_b_o=1./(0.000009./exp((v-26)/28.5)+0.014./(0.2+exp(-(v+70)/11)));