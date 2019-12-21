function tau_h_e=tau_h_e(v);
alpha_h=0.128*exp(-(v+50)/18);
beta_h=4./(1+exp(-(v+27)/5));
tau_h_e=1./(alpha_h+beta_h);
