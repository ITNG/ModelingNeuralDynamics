function beta_n_o=beta_n_o(v);
beta_n_o=0.0036*(35-v)./(1-exp(-(35-v)/12));