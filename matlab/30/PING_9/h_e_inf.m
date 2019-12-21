function h_e_inf=h_e_inf(v);
alpha_h=0.128*exp(-(v+50)/18);
beta_h=4./(1+exp(-(v+27)/5));
h_e_inf=alpha_h./(alpha_h+beta_h);