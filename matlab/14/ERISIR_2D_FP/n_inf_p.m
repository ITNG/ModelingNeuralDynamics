function n_inf_p=n_inf_p(v);
n_inf_p=(alpha_n(v)+beta_n(v)).*alpha_n_p(v);
n_inf_p=n_inf_p-alpha_n(v).*(alpha_n_p(v)+beta_n_p(v));
n_inf_p=n_inf_p./((alpha_n(v)+beta_n(v))).^2;