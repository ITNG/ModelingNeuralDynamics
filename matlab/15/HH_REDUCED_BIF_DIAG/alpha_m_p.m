function alpha_m_p=alpha_m_p(v);
num=(v+45)/10;
den=(1-exp(-(v+45)/10));
num_p=1/10;
den_p=exp(-(v+45)/10)/10;
alpha_m_p=(den.*num_p-num.*den_p)./(den.^2);


