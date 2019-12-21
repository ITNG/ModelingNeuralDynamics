function m_i_inf=m_i_inf(v);
alpha_m=0.1*(v+35)./(1-exp(-(v+35)/10));
beta_m=4*exp(-(v+60)/18);
m_i_inf=alpha_m./(alpha_m+beta_m);
