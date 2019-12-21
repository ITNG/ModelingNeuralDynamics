function tau_n_i=tau_n_i(v);
alpha_n=-0.01*(v+34)./(exp(-0.1*(v+34))-1);
beta_n=0.125*exp(-(v+44)/80);
tau_n_i=1./(alpha_n+beta_n);
phi=5; tau_n_i=tau_n_i/phi;
