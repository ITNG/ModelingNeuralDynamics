function beta_n=beta_n(v);

global phi;

beta_n=0.625*exp(-(v+44)/80);
beta_n=phi/5*beta_n;