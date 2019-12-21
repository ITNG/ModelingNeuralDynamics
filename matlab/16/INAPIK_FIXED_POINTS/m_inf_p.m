function m_inf_p=m_inf_p(v)
m_inf_p=-1./(1+exp((-20-v)/15))^2*exp((-20-v)/15)*(-1/15);