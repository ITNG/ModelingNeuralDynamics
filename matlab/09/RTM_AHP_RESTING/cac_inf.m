function cac_inf=cac_inf(v);
cac_inf=(120-v)./(1+exp(-(v+15)/5))*4/25;