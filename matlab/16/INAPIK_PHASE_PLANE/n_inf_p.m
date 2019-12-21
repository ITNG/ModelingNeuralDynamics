function n_inf_p=n_inf_p(v)
n_inf_p=-1./(1+exp((-25-v)/5)).^2.*exp((-25-v)/5)*(-1/5);