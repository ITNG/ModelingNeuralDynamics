function tau_h_i=tau_h_i(v)

phi=2.5;

alpha_h=0.07*exp(-(v+58)/20);
beta_h=1./(exp(-0.1*(v+28))+1);
tau_h_i=1./(alpha_h+beta_h); 
tau_h_i=tau_h_i/phi;
