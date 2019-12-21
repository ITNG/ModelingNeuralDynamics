function alpha_n_o=alpha_n_o(v)
alpha_n_o=0.018*(v-25)./(1-exp(-(v-25)/25));