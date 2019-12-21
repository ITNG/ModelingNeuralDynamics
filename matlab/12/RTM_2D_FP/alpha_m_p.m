function alpha_m_p=alpha_m_p(v);
num=0.32*(v+54);
den=1-exp(-(v+54)/4);
num_p=0.32;
den_p=exp(-(v+54)/4)/4;
alpha_m_p=(den.*num_p-num.*den_p)./(den.^2);


