function alpha_n=alpha_n(v)
alpha_n=0.018*(v-25)./(1-exp(-(v-25)/25));