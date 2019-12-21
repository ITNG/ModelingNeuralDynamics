function m_e_inf=m_e_inf(v);
alpha_m=0.32*(v+54)./(1-exp(-(v+54)/4));
beta_m=0.28*(v+27)./(exp((v+27)/5)-1);
m_e_inf=alpha_m./(alpha_m+beta_m);
