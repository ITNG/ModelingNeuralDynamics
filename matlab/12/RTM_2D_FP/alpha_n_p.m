function alpha_n_p=alpha_n_p(v);
num=0.032*(v+52); den=1-exp(-(v+52)/5);
num_p=0.032; den_p=exp(-(v+52)/5)/5;
alpha_n_p=(den.*num_p-num.*den_p)./(den.^2);
