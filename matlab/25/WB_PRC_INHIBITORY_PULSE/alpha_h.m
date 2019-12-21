function alpha_h=alpha_h(v); 

global phi;

alpha_h=0.35*exp(-(v+58)/20);
alpha_h=phi/5*alpha_h;
