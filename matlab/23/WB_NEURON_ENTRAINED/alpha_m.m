function alpha_m=alpha_m(v);
alpha_m=0.1*(v+35)./(1-exp(-(v+35)/10));

