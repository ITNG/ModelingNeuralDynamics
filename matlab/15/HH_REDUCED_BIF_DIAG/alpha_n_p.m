function alpha_n_p=alpha_n_p(v);
num=0.01*(-60.0d0-v); den=exp((-60-v)/10)-1;
num_p=-0.01; den_p=-(den+1)*0.1;
alpha_n_p=(den.*num_p-num.*den_p)./(den.^2);
