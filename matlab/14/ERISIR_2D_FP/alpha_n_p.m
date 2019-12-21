function alpha_n_p=alpha_n_p(v);
num=95-v; den=exp((95-v)/11.8)-1;
num_p=-1; den_p=-(den+1)/11.8;
alpha_n_p=(den.*num_p-num.*den_p)./(den.^2);

