function tau_r_modified=tau_r_modified(v);

tau_r=1./(exp(-14.59-0.086*v)+exp(-1.87+0.0701*v));
tau_r_modified=(((0.05*v+5).*(v<-80)+(v>=-80))).^2.*tau_r;


