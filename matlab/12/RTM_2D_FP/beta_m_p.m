function beta_m_p=beta_m_p(v);
num=0.28*(v+27);
den=exp((v+27)/5)-1;
num_p=0.28;
den_p=exp((v+27)/5)/5;
beta_m_p=(den.*num_p-num.*den_p)./(den.^2);

