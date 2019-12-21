function alpha_n=alpha_n(v);
alpha_n=0.032*(v+52)./(1-exp(-(v+52)/5));