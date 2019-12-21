function beta_n=beta_n(v);
beta_n=0.0036*(35-v)./(1-exp(-(35-v)/12));