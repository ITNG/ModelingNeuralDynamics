function h_i_inf=h_i_inf(v);
alpha_h=0.07*exp(-(v+58)/20);
beta_h=1./(exp(-0.1*(v+28))+1);
h_i_inf=alpha_h./(alpha_h+beta_h);