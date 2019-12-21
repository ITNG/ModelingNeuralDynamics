function beta_h=beta_h(v);

global phi;

beta_h=5./(exp(-0.1*(v+28))+1);
beta_h=phi/5*beta_h;
