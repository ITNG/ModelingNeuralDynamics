function g_0=g_0(phi);

phi_tilde=mod(phi,1);
g_0=phi_tilde.^2.*(1-phi_tilde);