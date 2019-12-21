function m_inf_p=m_inf_p(v);
m_inf_p=(alpha_m(v)+beta_m(v)).*alpha_m_p(v);
m_inf_p=m_inf_p-alpha_m(v).*(alpha_m_p(v)+beta_m_p(v));
m_inf_p=m_inf_p./((alpha_m(v)+beta_m(v))).^2;