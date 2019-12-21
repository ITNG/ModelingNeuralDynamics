function tau_b=tau_b(v)
tau_b=1./(0.000009./exp((v-26)/28.5)+0.014./(0.2+exp(-(v+70)/11)));