function alpha_n=alpha_n(v);

global phi;

alpha_n=0.05*(v+34)./(1-exp(-0.1*(v+34)));
alpha_n=phi/5*alpha_n;