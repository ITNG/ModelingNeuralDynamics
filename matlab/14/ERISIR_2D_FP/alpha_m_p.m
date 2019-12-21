function alpha_m_p=alpha_m_p(v);
num=40*(75.5-v);
den=exp((75.5-v)/13.5)-1;
num_p=-40;
den_p=-exp((75.5-v)/13.5)/13.5;
alpha_m_p=(den.*num_p-num.*den_p)./(den.^2);


